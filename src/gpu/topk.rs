//! Cached CUB workspace sizes and stream-ordered selection over existing distances.
use super::operations::allocate;
use anyhow::{ensure, Context, Result};
use cuda_async::{device_box::DeviceBox, device_operation::ExecutionContext, error::DeviceError};
use cuda_core::Context as CudaContext;
use std::{
    collections::HashMap,
    ffi::{c_char, c_int, c_void, CStr},
    mem, ptr,
};

// Handwritten C ABI: View matches C++'s iterator input, not a device allocation.
// Pointer-bearing views are constructed only for a query/submission, never cached.
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq)]
enum Layout {
    ReferenceRows = 0,
    QueryColumns = 1,
    DiagonalRows = 2,
    IntermediateRows = 3,
}
#[repr(C)]
struct View {
    keys: *const f32,
    ids: *const u32,
    width: usize,
    columns: usize,
    samples: usize,
    layout: Layout,
}
extern "C" {
    fn sketchlib_topk_limits(chunk: *mut usize, k: *mut usize);
    fn sketchlib_batched_topk(
        workspace: *mut c_void,
        bytes: *mut usize,
        input: *const View,
        chunks: usize,
        output_width: usize,
        k: c_int,
        segments: c_int,
        keys: *mut f32,
        ids: *mut u32,
        stream: *mut c_void,
    ) -> c_int;
    fn sketchlib_rows_topk(
        workspace: *mut c_void,
        bytes: *mut usize,
        input: *const View,
        rows: usize,
        k: usize,
        keys: *mut f32,
        ids: *mut u32,
        stream: *mut c_void,
    ) -> c_int;
    fn sketchlib_topk_error_string(status: c_int) -> *const c_char;
}
fn check(status: c_int) -> Result<(), DeviceError> {
    if status == 0 {
        return Ok(());
    }
    // CUDA runtime error strings are static and NUL terminated.
    let message = unsafe { CStr::from_ptr(sketchlib_topk_error_string(status)) }.to_string_lossy();
    Err(DeviceError::Internal(format!(
        "CUB tile selection: {message} (runtime error {status})"
    )))
}
#[derive(Clone, Copy)]
struct Limits {
    chunk: usize,
    k: usize,
}
impl Limits {
    fn compiled() -> Self {
        let mut limits = Self { chunk: 0, k: 0 };
        unsafe { sketchlib_topk_limits(&mut limits.chunk, &mut limits.k) };
        limits
    }
}

/// Workspace sizes are queried once per shape/device; cache hits copy scalars.
#[derive(Default)]
pub(super) struct SelectionCache {
    shapes: HashMap<(usize, usize, usize, bool), PreparedSelection>,
    limits: Option<Limits>,
}
impl SelectionCache {
    pub fn prepare(
        &mut self,
        context: &CudaContext,
        nr: usize,
        nq: usize,
        k: usize,
        diagonal: bool,
    ) -> Result<Shape> {
        let key = (nr, nq, k, diagonal);
        if let Some(selection) = self.shapes.get(&key) {
            return Ok(selection.shape);
        }
        let limits = *self.limits.get_or_insert_with(Limits::compiled);
        let selection = PreparedSelection::new(Shape::new(nr, nq, k, diagonal)?, limits, context)?;
        let shape = selection.shape;
        self.shapes.insert(key, selection);
        Ok(shape)
    }
    pub fn get(&self, nr: usize, nq: usize, k: usize, diagonal: bool) -> PreparedSelection {
        self.shapes[&(nr, nq, k, diagonal)]
    }
}
/// Candidate layout shared with C++ and the host merge: nr rows of ref_k pairs,
/// followed by nq rows of query_k pairs. Diagonals exclude self-pairs and have
/// no query section. IDs remain local to the opposite sketch block (the same
/// block for diagonals); candidates are unsorted and cutoff ties unspecified.
#[derive(Clone, Copy)]
pub(super) struct Shape {
    nr: usize,
    nq: usize,
    diagonal: bool,
    pub ref_k: usize,
    pub query_k: usize,
    pub count: usize,
    bytes: usize,
}

impl Shape {
    pub fn new(nr: usize, nq: usize, k: usize, diagonal: bool) -> Result<Self> {
        ensure!(
            nr > 0 && nq > 0 && k > 0,
            "CUB tile dimensions and k must be positive"
        );
        ensure!(
            nr <= u32::MAX as usize && nq <= u32::MAX as usize,
            "CUB block-local IDs exceed u32"
        );
        ensure!(!diagonal || nr == nq, "Diagonal CUB tile must be square");
        let ref_k = k.min(if diagonal { nr - 1 } else { nq });
        let query_k = if diagonal { 0 } else { k.min(nr) };
        let count = nr
            .checked_mul(ref_k)
            .and_then(|n| nq.checked_mul(query_k).and_then(|q| n.checked_add(q)))
            .context("CUB candidate count overflow")?;
        ensure!(
            count.checked_mul(mem::size_of::<f32>()).is_some(),
            "CUB candidate byte size overflow"
        );
        Ok(Self {
            nr,
            nq,
            diagonal,
            ref_k,
            query_k,
            count,
            bytes: 0,
        })
    }

    /// Device storage for candidate keys and IDs plus the queried CUB workspace.
    pub fn device_bytes(self) -> Result<usize> {
        self.count
            .checked_mul(mem::size_of::<f32>() + mem::size_of::<u32>())
            .and_then(|n| n.checked_add(self.bytes))
            .context("GPU selection byte size overflow")
    }
}

/// Initial row layout and scalar workspace requirements for one direction.
#[derive(Clone, Copy)]
struct Direction {
    layout: Layout,
    rows: usize,
    width: usize,
    k: usize,
    output_offset: usize,
    buffer_bytes: usize,
    scratch_bytes: usize,
}
fn reduced(width: usize, k: usize, chunk: usize) -> Result<usize> {
    (width / chunk)
        .checked_mul(k)
        .and_then(|n| n.checked_add((width % chunk).min(k)))
        .context("CUB reduction width overflow")
}
fn aligned(bytes: usize) -> Result<usize> {
    // CUDA allocation/scratch alignment, independent of CUB's selection bounds.
    const ALIGNMENT: usize = 256;
    bytes
        .checked_add(ALIGNMENT - 1)
        .map(|n| n & !(ALIGNMENT - 1))
        .context("CUB workspace alignment overflow")
}
impl Direction {
    fn new(shape: Shape, layout: Layout, limits: Limits) -> Result<Self> {
        let (rows, width, k, output_offset) = if layout == Layout::QueryColumns {
            // Shape already checked this product while counting candidates.
            (shape.nq, shape.nr, shape.query_k, shape.nr * shape.ref_k)
        } else {
            (
                shape.nr,
                if shape.diagonal {
                    shape.nr - 1
                } else {
                    shape.nq
                },
                shape.ref_k,
                0,
            )
        };
        let mut direction = Self {
            layout,
            rows,
            width,
            k,
            output_offset,
            buffer_bytes: 0,
            scratch_bytes: 0,
        };
        if k <= limits.k {
            let (_, output_width, _) = direction.batch_dimensions(width, limits)?;
            if width > limits.chunk {
                // The first reduction is largest; two key/ID pairs hold every pass.
                direction.buffer_bytes = aligned(
                    rows.checked_mul(output_width)
                        .and_then(|n| n.checked_mul(mem::size_of::<f32>()))
                        .context("CUB intermediate buffer size overflow")?,
                )?;
            }
        }
        direction.workspace_bytes()?;
        Ok(direction)
    }
    fn batch_dimensions(&self, width: usize, limits: Limits) -> Result<(usize, usize, c_int)> {
        let chunks = (width - 1) / limits.chunk + 1;
        let segments = self
            .rows
            .checked_mul(chunks)
            .and_then(|n| c_int::try_from(n).ok())
            .context("CUB segment count exceeds int")?;
        Ok((chunks, reduced(width, self.k, limits.chunk)?, segments))
    }
    fn view(&self, shape: Shape) -> View {
        View {
            keys: ptr::null(),
            ids: ptr::null(),
            width: self.width,
            columns: shape.nq,
            samples: shape.nr,
            layout: self.layout,
        }
    }
    fn workspace_bytes(&self) -> Result<usize> {
        self.buffer_bytes
            .checked_mul(4)
            .and_then(|n| n.checked_add(self.scratch_bytes.max(1)))
            .context("CUB workspace size overflow")
    }
    fn query(&mut self, shape: Shape, limits: Limits) -> Result<()> {
        let mut view = self.view(shape);
        if self.k > limits.k {
            check(unsafe {
                sketchlib_rows_topk(
                    ptr::null_mut(),
                    &mut self.scratch_bytes,
                    &view,
                    self.rows,
                    self.k,
                    ptr::null_mut(),
                    ptr::null_mut(),
                    ptr::null_mut(),
                )
            })?;
            return Ok(());
        }
        loop {
            let (chunks, output_width, segments) = self.batch_dimensions(view.width, limits)?;
            let mut bytes = 0;
            check(unsafe {
                sketchlib_batched_topk(
                    ptr::null_mut(),
                    &mut bytes,
                    &view,
                    chunks,
                    output_width,
                    self.k as c_int,
                    segments,
                    ptr::null_mut(),
                    ptr::null_mut(),
                    ptr::null_mut(),
                )
            })?;
            self.scratch_bytes = self.scratch_bytes.max(bytes);
            if view.width <= limits.chunk {
                break;
            }
            view.width = output_width;
            view.layout = Layout::IntermediateRows;
        }
        Ok(())
    }
}

/// Device-specific workspace sizes and initial layouts, without stored pass plans.
/// Copying this record needs no allocation or shared ownership of host metadata.
#[derive(Clone, Copy)]
pub(super) struct PreparedSelection {
    pub shape: Shape,
    limits: Limits,
    directions: [Option<Direction>; 2],
}
impl PreparedSelection {
    fn new(mut shape: Shape, limits: Limits, context: &CudaContext) -> Result<Self> {
        let mut directions = [None, None];
        if shape.ref_k != 0 {
            let layout = if shape.diagonal {
                Layout::DiagonalRows
            } else {
                Layout::ReferenceRows
            };
            directions[0] = Some(Direction::new(shape, layout, limits)?);
        }
        if shape.query_k != 0 {
            directions[1] = Some(Direction::new(shape, Layout::QueryColumns, limits)?);
        }
        // Queries may inspect device properties, even with null input pointers.
        context.bind_to_thread()?;
        shape.bytes = 1;
        for direction in directions.iter_mut().flatten() {
            direction.query(shape, limits)?;
            shape.bytes = shape.bytes.max(direction.workspace_bytes()?);
        }
        shape.device_bytes()?;
        Ok(Self {
            shape,
            limits,
            directions,
        })
    }
}

/// Enqueue cached selection on the tile stream, with three leased allocations.
///
/// # Safety
/// The execution context must be current, distances retained through completion,
/// and selection must have been prepared on this device for the distance layout.
pub(super) unsafe fn select(
    ctx: &ExecutionContext,
    distances: &DeviceBox<[f32]>,
    selection: PreparedSelection,
) -> Result<(DeviceBox<[f32]>, DeviceBox<[u32]>), DeviceError> {
    let device_keys = unsafe { allocate::<f32>(ctx, selection.shape.count)? };
    let device_ids = unsafe { allocate::<u32>(ctx, selection.shape.count)? };
    let workspace = unsafe { allocate::<u8>(ctx, selection.shape.bytes)? };
    let base = workspace.cu_deviceptr() as *mut u8;
    let stream = ctx.get_cuda_stream().cu_stream().cast();
    let limits = selection.limits;
    for direction in selection.directions.iter().flatten() {
        let keys = unsafe { (device_keys.cu_deviceptr() as *mut f32).add(direction.output_offset) };
        let ids = unsafe { (device_ids.cu_deviceptr() as *mut u32).add(direction.output_offset) };
        let mut view = direction.view(selection.shape);
        view.keys = distances.cu_deviceptr() as *const f32;
        if direction.k > limits.k {
            let mut bytes = direction.scratch_bytes;
            check(unsafe {
                sketchlib_rows_topk(
                    base.cast(),
                    &mut bytes,
                    &view,
                    direction.rows,
                    direction.k,
                    keys,
                    ids,
                    stream,
                )
            })?;
            continue;
        }
        let scratch = unsafe { base.add(4 * direction.buffer_bytes) }.cast();
        let mut slot = 0;
        loop {
            let (chunks, output_width, segments) =
                direction
                    .batch_dimensions(view.width, limits)
                    .map_err(|error| DeviceError::Internal(error.to_string()))?;
            let final_pass = view.width <= limits.chunk;
            let (out_keys, out_ids) = if final_pass {
                (keys, ids)
            } else {
                unsafe {
                    (
                        base.add(slot * 2 * direction.buffer_bytes).cast(),
                        base.add((slot * 2 + 1) * direction.buffer_bytes).cast(),
                    )
                }
            };
            let mut bytes = direction.scratch_bytes;
            check(unsafe {
                sketchlib_batched_topk(
                    scratch,
                    &mut bytes,
                    &view,
                    chunks,
                    output_width,
                    direction.k as c_int,
                    segments,
                    out_keys,
                    out_ids,
                    stream,
                )
            })?;
            if final_pass {
                break;
            }
            view.keys = out_keys;
            view.ids = out_ids;
            view.width = output_width;
            view.layout = Layout::IntermediateRows;
            slot ^= 1;
        }
    }
    Ok((device_keys, device_ids))
}

#[cfg(test)]
mod tests {
    use super::*;
    const LIMITS: Limits = Limits {
        chunk: 4096,
        k: 2048,
    };
    fn direction(width: usize, k: usize) -> Direction {
        Direction::new(
            Shape::new(3, width, k, false).unwrap(),
            Layout::ReferenceRows,
            LIMITS,
        )
        .unwrap()
    }

    #[test]
    fn reduction_dimensions_and_workspace_sizes() {
        let d = direction(20000, 2048);
        assert_eq!(d.batch_dimensions(20000, LIMITS).unwrap(), (5, 10240, 15));
        assert_eq!(d.batch_dimensions(10240, LIMITS).unwrap(), (3, 6144, 9));
        assert_eq!(d.batch_dimensions(6144, LIMITS).unwrap(), (2, 4096, 6));
        assert_eq!(d.batch_dimensions(4096, LIMITS).unwrap(), (1, 2048, 3));
        assert_eq!(d.buffer_bytes, 3 * 10240 * 4);
        assert_eq!(d.workspace_bytes().unwrap(), 4 * d.buffer_bytes + 1);
        for width in [4095, 4096] {
            let d = direction(width, 50);
            assert_eq!(d.batch_dimensions(width, LIMITS).unwrap(), (1, 50, 3));
            assert_eq!(d.buffer_bytes, 0);
        }
        assert!(direction(5003, 2048).buffer_bytes > 0);
        assert_eq!(direction(5003, 2049).buffer_bytes, 0);
    }
    #[test]
    fn shapes_clamp_each_direction_and_exclude_diagonal_self_pairs() {
        let s = Shape::new(3, 7, 50, false).unwrap();
        assert_eq!((s.ref_k, s.query_k, s.count), (7, 3, 42));
        let s = Shape::new(7, 7, 50, true).unwrap();
        assert_eq!((s.ref_k, s.query_k, s.count), (6, 0, 42));
        assert_eq!(Shape::new(1, 1, 50, true).unwrap().count, 0);
        assert!(Shape::new(3, 7, 2, true).is_err());
        assert!(Shape::new(0, 7, 2, false).is_err());
        assert!(Shape::new(3, 7, 0, false).is_err());
    }
    #[test]
    fn planning_rejects_overflows_before_workspace_queries() {
        assert!(aligned(usize::MAX).is_err());
        assert_eq!(aligned(257).unwrap(), 512);
        assert!(reduced(usize::MAX, usize::MAX, 4096).is_err());
        assert!(Direction::new(
            Shape::new(u32::MAX as usize, 8193, 2048, false).unwrap(),
            Layout::ReferenceRows,
            LIMITS
        )
        .is_err());
        assert!(Shape::new(
            u32::MAX as usize,
            u32::MAX as usize,
            u32::MAX as usize,
            false
        )
        .is_err());
    }
}
