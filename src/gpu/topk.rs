//! CUB workspace sizing and stream-ordered selection over existing distances.
use super::operations::allocate;
use anyhow::{ensure, Context, Result};
use cuda_async::{device_box::DeviceBox, device_operation::ExecutionContext, error::DeviceError};
use cuda_core::Context as CudaContext;
use std::{
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
}
#[repr(C)]
struct View {
    keys: *const f32,
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
pub(super) struct Limits {
    chunk: usize,
    k: usize,
}
impl Limits {
    pub fn compiled() -> Self {
        let mut limits = Self { chunk: 0, k: 0 };
        unsafe { sketchlib_topk_limits(&mut limits.chunk, &mut limits.k) };
        limits
    }
}

/// Candidate rows are packed by direction. IDs are local to the opposite block.
#[derive(Clone, Copy)]
pub(super) struct CandidateLayout {
    pub ref_stride: usize,
    pub query_stride: usize,
    pub count: usize,
}

#[derive(Clone, Copy)]
struct Batch {
    chunks: usize,
    segments: c_int,
}

#[derive(Clone, Copy)]
struct Direction {
    layout: Layout,
    rows: usize,
    width: usize,
    k: usize,
    stride: usize,
    batch: Option<Batch>,
}

impl Direction {
    fn new(layout: Layout, rows: usize, width: usize, k: usize, limits: Limits) -> Result<Self> {
        let k = k.min(width);
        let (stride, batch) = if k <= limits.k {
            let chunks = width.div_ceil(limits.chunk);
            let segments = rows
                .checked_mul(chunks)
                .and_then(|n| c_int::try_from(n).ok())
                .context("CUB segment count exceeds int")?;
            // Every full chunk contributes k winners. The last may be shorter.
            let stride = width / limits.chunk * k + (width % limits.chunk).min(k);
            (stride, Some(Batch { chunks, segments }))
        } else {
            (k, None)
        };
        Ok(Self {
            layout,
            rows,
            width,
            k,
            stride,
            batch,
        })
    }

    // A null workspace queries scratch size. Otherwise CUB writes this
    // direction's candidates into the supplied output section on the stream.
    unsafe fn dispatch(
        self,
        workspace: *mut c_void,
        bytes: &mut usize,
        view: &View,
        keys: *mut f32,
        ids: *mut u32,
        stream: *mut c_void,
    ) -> Result<(), DeviceError> {
        let status = unsafe {
            match self.batch {
                Some(batch) => sketchlib_batched_topk(
                    workspace,
                    bytes,
                    view,
                    batch.chunks,
                    self.stride,
                    self.k as c_int,
                    batch.segments,
                    keys,
                    ids,
                    stream,
                ),
                None => sketchlib_rows_topk(
                    workspace, bytes, view, self.rows, self.k, keys, ids, stream,
                ),
            }
        };
        check(status)
    }
}

/// Device-specific launches and shared scratch for the two selection directions.
#[derive(Clone, Copy)]
pub(super) struct PreparedSelection {
    pub layout: CandidateLayout,
    nr: usize,
    nq: usize,
    ref_rows: Option<Direction>,
    query_columns: Option<Direction>,
    scratch_bytes: usize,
}

impl PreparedSelection {
    fn plan(nr: usize, nq: usize, k: usize, diagonal: bool, limits: Limits) -> Result<Self> {
        ensure!(
            nr > 0 && nq > 0 && k > 0,
            "CUB tile dimensions and k must be positive"
        );
        ensure!(
            nr <= u32::MAX as usize && nq <= u32::MAX as usize,
            "CUB block-local IDs exceed u32"
        );
        ensure!(!diagonal || nr == nq, "Diagonal CUB tile must be square");

        let ref_rows = if diagonal {
            (nr > 1)
                .then(|| Direction::new(Layout::DiagonalRows, nr, nr - 1, k, limits))
                .transpose()?
        } else {
            Some(Direction::new(Layout::ReferenceRows, nr, nq, k, limits)?)
        };
        let query_columns = (!diagonal)
            .then(|| Direction::new(Layout::QueryColumns, nq, nr, k, limits))
            .transpose()?;
        let ref_stride = ref_rows.map_or(0, |d| d.stride);
        let query_stride = query_columns.map_or(0, |d| d.stride);
        let count = nr
            .checked_mul(ref_stride)
            .and_then(|n| nq.checked_mul(query_stride).and_then(|q| n.checked_add(q)))
            .context("CUB candidate count overflow")?;
        let selection = Self {
            layout: CandidateLayout {
                ref_stride,
                query_stride,
                count,
            },
            nr,
            nq,
            ref_rows,
            query_columns,
            scratch_bytes: 1,
        };
        selection.device_bytes()?;
        Ok(selection)
    }

    pub fn new(
        nr: usize,
        nq: usize,
        k: usize,
        diagonal: bool,
        limits: Limits,
        context: &CudaContext,
    ) -> Result<Self> {
        let mut selection = Self::plan(nr, nq, k, diagonal, limits)?;
        // CUB queries may inspect device properties, even with null inputs.
        context.bind_to_thread()?;
        for direction in [selection.ref_rows, selection.query_columns]
            .into_iter()
            .flatten()
        {
            let mut bytes = 0;
            unsafe {
                direction.dispatch(
                    ptr::null_mut(),
                    &mut bytes,
                    &selection.view(direction, ptr::null()),
                    ptr::null_mut(),
                    ptr::null_mut(),
                    ptr::null_mut(),
                )?;
            }
            selection.scratch_bytes = selection.scratch_bytes.max(bytes);
        }
        selection.device_bytes()?;
        Ok(selection)
    }

    pub fn device_bytes(self) -> Result<usize> {
        self.layout
            .count
            .checked_mul(mem::size_of::<f32>() + mem::size_of::<u32>())
            .and_then(|n| n.checked_add(self.scratch_bytes))
            .context("GPU selection byte size overflow")
    }

    fn view(self, direction: Direction, keys: *const f32) -> View {
        View {
            keys,
            width: direction.width,
            columns: self.nq,
            samples: self.nr,
            layout: direction.layout,
        }
    }
}

/// Enqueue prepared selection on the tile stream, with three leased allocations.
///
/// # Safety
/// The execution context must be current, distances retained through completion,
/// and selection must have been prepared on this device for the distance layout.
pub(super) unsafe fn select(
    ctx: &ExecutionContext,
    distances: &DeviceBox<[f32]>,
    selection: PreparedSelection,
) -> Result<(DeviceBox<[f32]>, DeviceBox<[u32]>), DeviceError> {
    let device_keys = unsafe { allocate::<f32>(ctx, selection.layout.count)? };
    let device_ids = unsafe { allocate::<u32>(ctx, selection.layout.count)? };
    let workspace = unsafe { allocate::<u8>(ctx, selection.scratch_bytes)? };
    let stream = ctx.get_cuda_stream().cu_stream().cast();
    let query_offset = selection.nr * selection.layout.ref_stride;
    for (direction, offset) in [
        (selection.ref_rows, 0),
        (selection.query_columns, query_offset),
    ] {
        if let Some(direction) = direction {
            let mut bytes = selection.scratch_bytes;
            unsafe {
                direction.dispatch(
                    workspace.cu_deviceptr() as *mut c_void,
                    &mut bytes,
                    &selection.view(direction, distances.cu_deviceptr() as *const f32),
                    (device_keys.cu_deviceptr() as *mut f32).add(offset),
                    (device_ids.cu_deviceptr() as *mut u32).add(offset),
                    stream,
                )?;
            }
        }
    }
    Ok((device_keys, device_ids))
}

#[cfg(test)]
mod tests {
    use super::*;
    const LIMITS: Limits = Limits {
        chunk: 8192,
        k: 2048,
    };

    #[test]
    fn candidate_strides_include_chunk_winners_and_row_fallback() {
        for (width, k, expected_stride, expected_segments) in [
            (8192, 50, 50, Some(3)),
            (8193, 50, 51, Some(6)),
            (10000, 50, 100, Some(6)),
            (25000, 50, 200, Some(12)),
            (8193, 2049, 2049, None),
        ] {
            let selection = PreparedSelection::plan(3, width, k, false, LIMITS).unwrap();
            let direction = selection.ref_rows.unwrap();
            assert_eq!(selection.layout.ref_stride, expected_stride);
            assert_eq!(selection.layout.query_stride, 3);
            assert_eq!(selection.layout.count, 3 * expected_stride + width * 3);
            assert_eq!(direction.batch.map(|b| b.segments), expected_segments);
        }
    }

    #[test]
    fn selection_clamps_k_and_excludes_diagonal_self_pairs() {
        let layout = PreparedSelection::plan(7, 7, 50, true, LIMITS)
            .unwrap()
            .layout;
        assert_eq!(
            (layout.ref_stride, layout.query_stride, layout.count),
            (6, 0, 42)
        );
        assert_eq!(
            PreparedSelection::plan(1, 1, 50, true, LIMITS)
                .unwrap()
                .layout
                .count,
            0
        );
        assert!(PreparedSelection::plan(3, 7, 2, true, LIMITS).is_err());
        assert!(PreparedSelection::plan(0, 7, 2, false, LIMITS).is_err());
        assert!(PreparedSelection::plan(3, 7, 0, false, LIMITS).is_err());
    }

    #[test]
    fn planning_rejects_overflows_before_workspace_queries() {
        assert!(PreparedSelection::plan(u32::MAX as usize, 8193, 2048, false, LIMITS).is_err());
        assert!(PreparedSelection::plan(
            u32::MAX as usize,
            u32::MAX as usize,
            u32::MAX as usize,
            false,
            LIMITS,
        )
        .is_err());
    }
}
