//! Private CUB dispatch and device candidate storage.
use anyhow::{ensure, Context, Result};
use cuda_async::device_operation::BoxedDeviceOp;
use cuda_async::{
    device_box::DeviceBox,
    device_operation::{DeviceOperation, ExecutionContext},
    error::DeviceError,
};
use cuda_core::memory::{malloc_async, memcpy_dtoh_async};
use cuda_core::{Context as CudaContext, DriverError, PinnedHostBuffer};
use std::future::IntoFuture;
use std::{
    ffi::{c_char, c_int, c_void, CStr},
    mem, ptr,
    sync::Arc,
};

extern "C" {
    fn sketchlib_tile_topk(
        workspace: *mut c_void,
        bytes: *mut usize,
        distances: *const f32,
        keys: *mut f32,
        ids: *mut u32,
        nr: usize,
        nq: usize,
        k: usize,
        diagonal: c_int,
        stream: *mut c_void,
    ) -> c_int;
    fn sketchlib_init_candidates(
        keys: *mut f32,
        ids: *mut u32,
        n: usize,
        block: usize,
        k: usize,
        width: usize,
        stream: *mut c_void,
    ) -> c_int;
    fn sketchlib_scatter_candidates(
        keys: *const f32,
        ids: *const u32,
        all_keys: *mut f32,
        all_ids: *mut u32,
        nr: usize,
        nq: usize,
        k: usize,
        ref_start: usize,
        query_start: usize,
        block: usize,
        width: usize,
        stream: *mut c_void,
    ) -> c_int;
    fn sketchlib_global_topk(
        workspace: *mut c_void,
        bytes: *mut usize,
        keys: *const f32,
        ids: *const u32,
        out_keys: *mut f32,
        out_ids: *mut u32,
        rows: usize,
        width: usize,
        k: usize,
        stream: *mut c_void,
    ) -> c_int;
    fn sketchlib_topk_error_string(status: c_int) -> *const c_char;
}

fn check(status: c_int) -> Result<(), DeviceError> {
    if status == 0 {
        return Ok(());
    }
    // CUDA runtime error strings are static, NUL-terminated strings.
    let message = unsafe { CStr::from_ptr(sketchlib_topk_error_string(status)) }.to_string_lossy();
    Err(DeviceError::Internal(format!(
        "CUB selection: {message} (runtime error {status})"
    )))
}

#[derive(Clone, Copy)]
pub(super) struct Shape {
    nr: usize,
    nq: usize,
    k: usize,
    diagonal: bool,
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
            k,
            diagonal,
            count,
            bytes: 0,
        })
    }

    pub fn query(mut self) -> Result<Self> {
        // Query mode never dereferences these device pointers or enqueues work.
        check(unsafe {
            sketchlib_tile_topk(
                ptr::null_mut(),
                &mut self.bytes,
                ptr::null(),
                ptr::null_mut(),
                ptr::null_mut(),
                self.nr,
                self.nq,
                self.k,
                self.diagonal.into(),
                ptr::null_mut(),
            )
        })?;
        Ok(self)
    }
}

/// Every row has a separate segment for each neighbour block. Only the own-
/// block segment can have padding (one entry when that block has <= k samples).
#[derive(Clone, Copy, Debug)]
pub(super) struct Layout {
    n: usize,
    block: usize,
    k: usize,
    width: usize,
    count: usize,
    pub bytes: usize,
    output_count: usize,
}

impl Layout {
    pub fn new(n: usize, block: usize, k: usize) -> Result<Self> {
        ensure!(
            n > 0 && block > 0 && k > 0 && k < n,
            "Invalid GPU kNN layout"
        );
        ensure!(
            n <= u32::MAX as usize,
            "GPU global neighbour IDs exceed u32"
        );
        let blocks = n.div_ceil(block);
        let width = (blocks - 1)
            .checked_mul(k.min(block))
            .and_then(|prefix| prefix.checked_add(k.min(n - (blocks - 1) * block)))
            .context("GPU candidate row width overflow")?;
        let count = n
            .checked_mul(width)
            .context("GPU candidate count overflow")?;
        let bytes = count
            .checked_mul(8)
            .context("GPU candidate byte size overflow")?;
        let output_count = n
            .checked_mul(k)
            .context("GPU final candidate count overflow")?;
        ensure!(
            output_count
                .checked_mul(mem::size_of::<
                    crate::distances::distance_matrix::SparseJaccard,
                >())
                .is_some(),
            "GPU sparse output byte size overflow"
        );
        Ok(Self {
            n,
            block,
            k,
            width,
            count,
            bytes,
            output_count,
        })
    }
}

pub(super) struct Candidates {
    layout: Layout,
    keys: DeviceBox<[f32]>,
    ids: DeviceBox<[u32]>,
}

pub(super) enum TileOutput {
    Dense(Arc<PinnedHostBuffer<f32>>),
    Candidates(Arc<Candidates>),
}

impl TileOutput {
    pub fn new(
        ctx: &Arc<CudaContext>,
        max_distances: usize,
        candidates: Option<&Arc<Candidates>>,
    ) -> Result<Self, DriverError> {
        Ok(match candidates {
            None => Self::Dense(Arc::new(PinnedHostBuffer::zeroed(ctx, max_distances)?)),
            Some(candidates) => Self::Candidates(Arc::clone(candidates)),
        })
    }
    pub fn dense(&self) -> &[f32] {
        match self {
            Self::Dense(data) => data,
            _ => unreachable!("dense consumer received candidates"),
        }
    }
}

/// Await once before dispatching tiles, establishing allocation/initialisation
/// completion before any other stream uses the arrays.
pub(super) struct Prepare(pub Layout);
impl DeviceOperation for Prepare {
    type Output = Arc<Candidates>;
    unsafe fn execute(
        self,
        ctx: &ExecutionContext,
    ) -> Result<<Self as DeviceOperation>::Output, DeviceError> {
        let layout = self.0;
        let keys = unsafe { allocate::<f32>(ctx, layout.count)? };
        let ids = unsafe { allocate::<u32>(ctx, layout.count)? };
        check(unsafe {
            sketchlib_init_candidates(
                keys.cu_deviceptr() as *mut f32,
                ids.cu_deviceptr() as *mut u32,
                layout.n,
                layout.block,
                layout.k,
                layout.width,
                ctx.get_cuda_stream().cu_stream().cast(),
            )
        })?;
        Ok(Arc::new(Candidates { layout, keys, ids }))
    }
}
impl IntoFuture for Prepare {
    type Output = Result<Arc<Candidates>, DeviceError>;
    type IntoFuture = <BoxedDeviceOp<Arc<Candidates>> as IntoFuture>::IntoFuture;
    fn into_future(self) -> Self::IntoFuture {
        self.boxed().into_future()
    }
}

/// This operation runs only after all tile futures have completed successfully.
/// It owns pinned destinations and retains all device storage through download.
pub(super) struct Select {
    candidates: Arc<Candidates>,
    keys: Arc<PinnedHostBuffer<f32>>,
    ids: Arc<PinnedHostBuffer<u32>>,
    bytes: usize,
}
pub(super) type Selected = (Arc<PinnedHostBuffer<f32>>, Arc<PinnedHostBuffer<u32>>);
impl Select {
    pub fn new(candidates: Arc<Candidates>, ctx: &Arc<CudaContext>) -> Result<Self> {
        let l = candidates.layout;
        let mut bytes = 0;
        check(unsafe {
            sketchlib_global_topk(
                ptr::null_mut(),
                &mut bytes,
                ptr::null(),
                ptr::null(),
                ptr::null_mut(),
                ptr::null_mut(),
                l.n,
                l.width,
                l.k,
                ptr::null_mut(),
            )
        })?;
        log::info!(
            "Final GPU kNN selection: {} scratch bytes, {} result bytes",
            bytes,
            l.output_count * 8
        );
        Ok(Self {
            candidates,
            keys: Arc::new(PinnedHostBuffer::zeroed(ctx, l.output_count)?),
            ids: Arc::new(PinnedHostBuffer::zeroed(ctx, l.output_count)?),
            bytes,
        })
    }
}
impl DeviceOperation for Select {
    type Output = Selected;
    unsafe fn execute(
        mut self,
        ctx: &ExecutionContext,
    ) -> Result<<Self as DeviceOperation>::Output, DeviceError> {
        let l = self.candidates.layout;
        let stream = ctx.get_cuda_stream().cu_stream();
        ctx.retain(Arc::clone(&self.candidates))?;
        let host_keys = destination(ctx, &mut self.keys)?;
        let host_ids = destination(ctx, &mut self.ids)?;
        let keys = unsafe { allocate::<f32>(ctx, l.output_count)? };
        let ids = unsafe { allocate::<u32>(ctx, l.output_count)? };
        let workspace = unsafe { allocate::<u8>(ctx, self.bytes)? };
        check(unsafe {
            sketchlib_global_topk(
                workspace.cu_deviceptr() as *mut c_void,
                &mut self.bytes,
                self.candidates.keys.cu_deviceptr() as *const f32,
                self.candidates.ids.cu_deviceptr() as *const u32,
                keys.cu_deviceptr() as *mut f32,
                ids.cu_deviceptr() as *mut u32,
                l.n,
                l.width,
                l.k,
                stream.cast(),
            )
        })?;
        unsafe {
            memcpy_dtoh_async(host_keys, keys.cu_deviceptr(), l.output_count * 4, stream)?;
            memcpy_dtoh_async(host_ids, ids.cu_deviceptr(), l.output_count * 4, stream)?;
        }
        Ok((self.keys, self.ids))
    }
}
impl IntoFuture for Select {
    type Output = Result<Selected, DeviceError>;
    type IntoFuture = <BoxedDeviceOp<Selected> as IntoFuture>::IntoFuture;
    fn into_future(self) -> Self::IntoFuture {
        self.boxed().into_future()
    }
}

/// Register an owner before enqueueing raw asynchronous accesses to it.
pub(super) struct KeepAlive<T: Send + 'static>(pub T);
impl<T: Send + 'static> DeviceOperation for KeepAlive<T> {
    type Output = ();
    unsafe fn execute(self, ctx: &ExecutionContext) -> Result<(), DeviceError> {
        ctx.retain(self.0)
    }
}

impl<T: Send + 'static> IntoFuture for KeepAlive<T> {
    type Output = Result<(), DeviceError>;
    type IntoFuture = <BoxedDeviceOp<()> as IntoFuture>::IntoFuture;
    fn into_future(self) -> Self::IntoFuture {
        self.boxed().into_future()
    }
}

pub(super) struct Finish {
    pub distances: DeviceBox<[f32]>,
    pub target: TileOutput,
    pub selection: Option<Shape>,
    pub sample_offsets: (usize, usize),
}

// Allocate only on this operation's stream; leases protect raw CUB accesses.
unsafe fn allocate<T: cuda_core::DeviceCopy + Send>(
    ctx: &ExecutionContext,
    n: usize,
) -> Result<DeviceBox<[T]>, DeviceError> {
    let bytes = n
        .checked_mul(mem::size_of::<T>())
        .ok_or_else(|| DeviceError::Internal("CUB allocation overflow".into()))?;
    let stream = ctx.get_cuda_stream();
    let ptr = unsafe { malloc_async(stream.cu_stream(), bytes) }.map_err(|error| {
        DeviceError::Internal(format!("CUB allocation of {bytes} bytes failed: {error}"))
    })?;
    let buffer = unsafe { DeviceBox::from_raw_parts(ptr, n, Arc::clone(stream)) };
    ctx.retain_buffer(buffer.storage().clone(), true)?;
    Ok(buffer)
}

/// Obtain the DMA pointer exclusively before retaining a second owner. No Rust
/// reads or writes occur until the completion lease releases that owner.
fn destination<T: cuda_core::DeviceCopy + Send + Sync + 'static>(
    ctx: &ExecutionContext,
    buffer: &mut Arc<PinnedHostBuffer<T>>,
) -> Result<*mut T, DeviceError> {
    let pointer = Arc::get_mut(buffer)
        .ok_or_else(|| DeviceError::Internal("Pinned output is still in use".into()))?
        .as_mut_ptr();
    ctx.retain(Arc::clone(buffer))?;
    Ok(pointer)
}

impl DeviceOperation for Finish {
    type Output = (TileOutput, DeviceBox<[f32]>);
    unsafe fn execute(
        mut self,
        ctx: &ExecutionContext,
    ) -> Result<(TileOutput, DeviceBox<[f32]>), DeviceError> {
        let stream = ctx.get_cuda_stream().cu_stream();
        ctx.retain_buffer(self.distances.storage().clone(), false)?;
        match (&mut self.target, self.selection) {
            (TileOutput::Dense(host), None) => {
                let pointer = destination(ctx, host)?;
                unsafe {
                    memcpy_dtoh_async(
                        pointer,
                        self.distances.cu_deviceptr(),
                        self.distances.len() * mem::size_of::<f32>(),
                        stream,
                    )?;
                }
            }
            (TileOutput::Candidates(candidates), Some(mut shape)) => {
                // The layout gives each tile disjoint destinations. Retain the
                // owner without claiming a whole-allocation write lease, which
                // would serialize unrelated tiles. No reads occur until every
                // tile's completion event has been awaited.
                ctx.retain(Arc::clone(candidates))?;
                let device_keys = unsafe { allocate::<f32>(ctx, shape.count)? };
                let device_ids = unsafe { allocate::<u32>(ctx, shape.count)? };
                let workspace = unsafe { allocate::<u8>(ctx, shape.bytes)? };
                check(unsafe {
                    sketchlib_tile_topk(
                        workspace.cu_deviceptr() as *mut c_void,
                        &mut shape.bytes,
                        self.distances.cu_deviceptr() as *const f32,
                        device_keys.cu_deviceptr() as *mut f32,
                        device_ids.cu_deviceptr() as *mut u32,
                        shape.nr,
                        shape.nq,
                        shape.k,
                        shape.diagonal.into(),
                        stream.cast(),
                    )
                })?;
                let l = candidates.layout;
                check(unsafe {
                    sketchlib_scatter_candidates(
                        device_keys.cu_deviceptr() as *const f32,
                        device_ids.cu_deviceptr() as *const u32,
                        candidates.keys.cu_deviceptr() as *mut f32,
                        candidates.ids.cu_deviceptr() as *mut u32,
                        shape.nr,
                        shape.nq,
                        shape.k,
                        self.sample_offsets.0,
                        self.sample_offsets.1,
                        l.block,
                        l.width,
                        stream.cast(),
                    )
                })?;
            }
            _ => return Err(DeviceError::Internal("Tile output mode mismatch".into())),
        }
        Ok((self.target, self.distances))
    }
}

impl IntoFuture for Finish {
    type Output = Result<(TileOutput, DeviceBox<[f32]>), DeviceError>;
    type IntoFuture = <BoxedDeviceOp<(TileOutput, DeviceBox<[f32]>)> as IntoFuture>::IntoFuture;
    fn into_future(self) -> Self::IntoFuture {
        self.boxed().into_future()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn candidate_layout_caps_short_blocks() {
        // Blocks of 3, 3, 1 reserve 2, 2, 1 entries in each sample's row.
        let layout = Layout::new(7, 3, 2).unwrap();
        assert_eq!((layout.width, layout.count, layout.bytes), (5, 35, 280));
        assert_eq!(layout.output_count, 14);
        // k larger than every block never reserves more than n entries per row.
        let layout = Layout::new(7, 3, 6).unwrap();
        assert_eq!((layout.width, layout.count), (7, 49));
        assert_eq!(Layout::new(7, 1, 6).unwrap().width, 7);
        assert_eq!(Layout::new(7, 100, 6).unwrap().width, 6);
    }

    #[test]
    fn candidate_layout_accounts_for_diagonal_padding() {
        for (n, block, k) in [(7, 3, 2), (7, 3, 6), (7, 1, 6), (7, 100, 6)] {
            let layout = Layout::new(n, block, k).unwrap();
            for owner in (0..n).step_by(block) {
                let own_size = block.min(n - owner);
                let actual: usize = (0..n)
                    .step_by(block)
                    .map(|start| k.min(block.min(n - start) - usize::from(start == owner)))
                    .sum();
                assert_eq!(layout.width - actual, usize::from(own_size <= k));
            }
        }
    }

    #[test]
    fn candidate_layout_rejects_invalid_sizes_and_overflow() {
        for (n, block, k) in [(0, 1, 1), (7, 0, 1), (7, 3, 0), (7, 3, 7)] {
            assert!(Layout::new(n, block, k).is_err());
        }
        assert!(Layout::new(u32::MAX as usize + 1, 1, 1).is_err());
        assert!(Layout::new(u32::MAX as usize, 1, u32::MAX as usize - 1).is_err());
    }
}
