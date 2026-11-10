//! Private CUB host dispatch and stream-owned output buffers.
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
    fn sketchlib_topk_error_string(status: c_int) -> *const c_char;
}

fn check(status: c_int) -> Result<(), DeviceError> {
    if status == 0 {
        return Ok(());
    }
    // CUDA runtime error strings are static, NUL-terminated strings.
    let message = unsafe { CStr::from_ptr(sketchlib_topk_error_string(status)) }.to_string_lossy();
    Err(DeviceError::Internal(format!(
        "CUB tile selection: {message} (runtime error {status})"
    )))
}

#[derive(Clone, Copy)]
pub(super) struct Shape {
    nr: usize,
    nq: usize,
    k: usize,
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
            k,
            diagonal,
            ref_k,
            query_k,
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

pub(super) enum HostOutput {
    Dense(Arc<PinnedHostBuffer<f32>>),
    Candidates {
        keys: Arc<PinnedHostBuffer<f32>>,
        ids: Arc<PinnedHostBuffer<u32>>,
    },
}

impl HostOutput {
    pub fn new(
        ctx: &Arc<CudaContext>,
        max_distances: usize,
        candidates: Option<usize>,
    ) -> Result<Self, DriverError> {
        Ok(match candidates {
            None => Self::Dense(Arc::new(PinnedHostBuffer::zeroed(ctx, max_distances)?)),
            Some(n) => Self::Candidates {
                keys: Arc::new(PinnedHostBuffer::zeroed(ctx, n)?),
                ids: Arc::new(PinnedHostBuffer::zeroed(ctx, n)?),
            },
        })
    }
    pub fn dense(&self) -> &[f32] {
        match self {
            Self::Dense(data) => data,
            _ => unreachable!("dense consumer received candidates"),
        }
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
    pub host: HostOutput,
    pub selection: Option<Shape>,
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
    let ptr = unsafe { malloc_async(stream.cu_stream(), bytes)? };
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
    type Output = (HostOutput, DeviceBox<[f32]>);
    unsafe fn execute(
        mut self,
        ctx: &ExecutionContext,
    ) -> Result<(HostOutput, DeviceBox<[f32]>), DeviceError> {
        let stream = ctx.get_cuda_stream().cu_stream();
        ctx.retain_buffer(self.distances.storage().clone(), false)?;
        match (&mut self.host, self.selection) {
            (HostOutput::Dense(host), None) => {
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
            (HostOutput::Candidates { keys, ids }, Some(mut shape)) => {
                let host_keys = destination(ctx, keys)?;
                let host_ids = destination(ctx, ids)?;
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
                unsafe {
                    memcpy_dtoh_async(
                        host_keys,
                        device_keys.cu_deviceptr(),
                        shape.count * mem::size_of::<f32>(),
                        stream,
                    )?;
                    memcpy_dtoh_async(
                        host_ids,
                        device_ids.cu_deviceptr(),
                        shape.count * mem::size_of::<u32>(),
                        stream,
                    )?;
                }
                // Leases keep these allocations alive until the entire stream
                // submission, including both D2H copies, has completed.
            }
            _ => return Err(DeviceError::Internal("Tile output mode mismatch".into())),
        }
        Ok((self.host, self.distances))
    }
}

impl IntoFuture for Finish {
    type Output = Result<(HostOutput, DeviceBox<[f32]>), DeviceError>;
    type IntoFuture = <BoxedDeviceOp<(HostOutput, DeviceBox<[f32]>)> as IntoFuture>::IntoFuture;
    fn into_future(self) -> Self::IntoFuture {
        self.boxed().into_future()
    }
}
