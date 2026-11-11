//! Pinned tile output and stream-ordered transfer to host consumers.
use super::topk;
use cuda_async::{
    device_box::DeviceBox,
    device_operation::{BoxedDeviceOp, DeviceOperation, ExecutionContext},
    error::DeviceError,
};
use cuda_core::{memory::memcpy_dtoh_async, Context as CudaContext, DriverError, PinnedHostBuffer};
use std::{future::IntoFuture, mem, sync::Arc};

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
    pub selection: Option<topk::Shape>,
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
            (HostOutput::Candidates { keys, ids }, Some(shape)) => {
                let host_keys = destination(ctx, keys)?;
                let host_ids = destination(ctx, ids)?;
                let (device_keys, device_ids) =
                    unsafe { topk::select(ctx, &self.distances, shape)? };
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
    use crate::gpu::topk::SelectionCache;
    use cuda_async::device_context::{global_policy, init_device_contexts, with_cuda_context};
    use std::ptr;

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
