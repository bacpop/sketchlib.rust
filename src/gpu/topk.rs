//! Private CUB host dispatch and stream-owned output buffers.
use super::operations::allocate;
use anyhow::{ensure, Context, Result};
use cuda_async::device_operation::BoxedDeviceOp;
use cuda_async::{
    device_box::DeviceBox,
    device_operation::{DeviceOperation, ExecutionContext},
    error::DeviceError,
};
use cuda_core::memory::memcpy_dtoh_async;
use cuda_core::{Context as CudaContext, DriverError, PinnedHostBuffer};
use std::future::IntoFuture;
use std::{
    collections::HashMap,
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

/// Prepared selection shapes for one device, reused during planning and dispatch.
#[derive(Default)]
pub(super) struct SelectionCache {
    shapes: HashMap<(usize, usize, usize, bool), Shape>,
}

impl SelectionCache {
    /// Query on a cache miss. Each cache must be used with only one device.
    pub fn prepare(
        &mut self,
        context: &CudaContext,
        nr: usize,
        nq: usize,
        k: usize,
        diagonal: bool,
    ) -> Result<Shape> {
        let key = (nr, nq, k, diagonal);
        if let Some(shape) = self.shapes.get(&key) {
            return Ok(*shape);
        }
        let shape = Shape::new(nr, nq, k, diagonal)?.query_on(context)?;
        self.shapes.insert(key, shape);
        Ok(shape)
    }

    /// Retrieve a shape prepared during memory planning; no CUDA calls occur here.
    pub fn get(&self, nr: usize, nq: usize, k: usize, diagonal: bool) -> Shape {
        self.shapes[&(nr, nq, k, diagonal)]
    }
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

    /// Device storage for candidate keys and IDs plus the queried CUB workspace.
    pub fn device_bytes(self) -> Result<usize> {
        self.count
            .checked_mul(mem::size_of::<f32>() + mem::size_of::<u32>())
            .and_then(|n| n.checked_add(self.bytes))
            .context("GPU selection byte size overflow")
    }

    pub fn query_on(mut self, context: &CudaContext) -> Result<Self> {
        context.bind_to_thread()?;
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

pub(super) struct Finish {
    pub distances: DeviceBox<[f32]>,
    pub host: HostOutput,
    pub selection: Option<Shape>,
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
        // Runtime CUB uses the thread's current Driver API context, not just
        // the stream argument. Resource cleanup may have changed that context.
        ctx.context().bind_to_thread()?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use cuda_async::device_context::{global_policy, init_device_contexts, with_cuda_context};

    #[test]
    #[ignore = "requires three CUDA GPUs; run explicitly for multi-GPU validation"]
    fn cub_and_raw_operations_bind_their_own_context() {
        use cuda_core::IntoResult;
        assert!(CudaContext::device_count().unwrap() >= 3);
        init_device_contexts(0, 3).unwrap();
        let contexts: Vec<_> = (0..3)
            .map(|id| with_cuda_context(id, Arc::clone).unwrap())
            .collect();
        let data: Vec<_> = (0..35).map(|i| i as f32 / 100.0).collect();
        let mut jobs = Vec::new();
        for (id, context) in contexts.iter().enumerate() {
            let foreign = &contexts[(id + 1) % 3];
            let foreign_buffer = PinnedHostBuffer::<u32>::zeroed(foreign, 1).unwrap();
            context.bind_to_thread().unwrap();
            drop(foreign_buffer); // Deliberately changes the current context.
            let mut selections = SelectionCache::default();
            selections.prepare(context, 5, 7, 3, false).unwrap();
            let shape = selections.get(5, 7, 3, false);
            let mut current = ptr::null_mut();
            unsafe {
                cuda_core::sys::cuCtxGetCurrent(&mut current)
                    .result()
                    .unwrap();
            }
            assert_eq!(current, context.cu_ctx());
            let policy = global_policy(id).unwrap();
            let stream = policy.next_stream().unwrap();
            foreign.bind_to_thread().unwrap();
            let distances = crate::gpu::operations::h2d(&data).sync_on(&stream).unwrap();
            let host = HostOutput::new(context, 0, Some(shape.count)).unwrap();
            let job = Finish {
                distances,
                host,
                selection: Some(shape),
            }
            .schedule(policy.as_ref())
            .unwrap();
            foreign.bind_to_thread().unwrap();
            jobs.push(job);
        }
        crate::gpu::run_gpu_batches(async {
            let mut pending = tokio::task::JoinSet::new();
            for job in jobs {
                pending.spawn(job);
            }
            while let Some(result) = pending.join_next().await {
                let (host, _distances) = result.unwrap().unwrap();
                let HostOutput::Candidates { keys, ids } = host else {
                    panic!("Expected candidates")
                };
                for row in 0..5 {
                    let mut chosen = Vec::new();
                    for j in 0..3 {
                        let offset = row * 3 + j;
                        let id = ids[offset] as usize;
                        assert_eq!(keys[offset], data[row * 7 + id]);
                        chosen.push(id);
                    }
                    chosen.sort_unstable();
                    assert_eq!(chosen, [0, 1, 2]);
                }
                for col in 0..7 {
                    let mut chosen = Vec::new();
                    for j in 0..3 {
                        let offset = 15 + col * 3 + j;
                        let id = ids[offset] as usize;
                        assert_eq!(keys[offset], data[id * 7 + col]);
                        chosen.push(id);
                    }
                    chosen.sort_unstable();
                    assert_eq!(chosen, [0, 1, 2]);
                }
            }
            Ok(())
        })
        .unwrap();
    }
}
