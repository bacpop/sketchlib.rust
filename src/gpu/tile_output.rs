//! Pinned tile output and stream-ordered transfer to host consumers.
use super::topk;
use cuda_async::{
    device_box::DeviceBox,
    device_operation::{BoxedDeviceOp, DeviceOperation, ExecutionContext},
    error::DeviceError,
};
use cuda_core::{memory::memcpy_dtoh_async, Context as CudaContext, DriverError, PinnedHostBuffer};
use std::{future::IntoFuture, mem, sync::Arc};

/// Pinned buffers, with a candidate layout only after a tile fills them.
pub(super) enum TileOutput<L> {
    Dense(Arc<PinnedHostBuffer<f32>>),
    Candidates(CandidateOutput<L>),
}

pub(super) struct CandidateOutput<L> {
    pub keys: Arc<PinnedHostBuffer<f32>>,
    pub ids: Arc<PinnedHostBuffer<u32>>,
    pub layout: L,
}

pub(super) type HostOutput = TileOutput<()>;
pub(super) type CompletedOutput = TileOutput<topk::CandidateLayout>;
pub(super) type CompletedCandidates = CandidateOutput<topk::CandidateLayout>;

impl TileOutput<()> {
    pub fn new(
        ctx: &Arc<CudaContext>,
        max_distances: usize,
        candidates: Option<usize>,
    ) -> Result<Self, DriverError> {
        Ok(match candidates {
            None => Self::Dense(Arc::new(PinnedHostBuffer::zeroed(ctx, max_distances)?)),
            Some(n) => Self::Candidates(CandidateOutput {
                keys: Arc::new(PinnedHostBuffer::zeroed(ctx, n)?),
                ids: Arc::new(PinnedHostBuffer::zeroed(ctx, n)?),
                layout: (),
            }),
        })
    }
}

impl CandidateOutput<topk::CandidateLayout> {
    pub fn into_buffers(self) -> HostOutput {
        HostOutput::Candidates(CandidateOutput {
            keys: self.keys,
            ids: self.ids,
            layout: (),
        })
    }
}

pub(super) struct Finish {
    pub distances: DeviceBox<[f32]>,
    pub host: HostOutput,
    pub selection: Option<topk::PreparedSelection>,
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
    type Output = (CompletedOutput, DeviceBox<[f32]>);
    unsafe fn execute(
        self,
        ctx: &ExecutionContext,
    ) -> Result<(CompletedOutput, DeviceBox<[f32]>), DeviceError> {
        // Runtime CUB uses the thread's current Driver API context, not just
        // the stream argument. Resource cleanup may have changed that context.
        ctx.context().bind_to_thread()?;
        let stream = ctx.get_cuda_stream().cu_stream();
        ctx.retain_buffer(self.distances.storage().clone(), false)?;
        let output = match (self.host, self.selection) {
            (HostOutput::Dense(mut host), None) => {
                let pointer = destination(ctx, &mut host)?;
                unsafe {
                    memcpy_dtoh_async(
                        pointer,
                        self.distances.cu_deviceptr(),
                        self.distances.len() * mem::size_of::<f32>(),
                        stream,
                    )?;
                }
                CompletedOutput::Dense(host)
            }
            (
                HostOutput::Candidates(CandidateOutput {
                    mut keys, mut ids, ..
                }),
                Some(selection),
            ) => {
                let count = selection.layout.count;
                let host_keys = destination(ctx, &mut keys)?;
                let host_ids = destination(ctx, &mut ids)?;
                let (device_keys, device_ids) =
                    unsafe { topk::select(ctx, &self.distances, selection)? };
                unsafe {
                    memcpy_dtoh_async(
                        host_keys,
                        device_keys.cu_deviceptr(),
                        count * mem::size_of::<f32>(),
                        stream,
                    )?;
                    memcpy_dtoh_async(
                        host_ids,
                        device_ids.cu_deviceptr(),
                        count * mem::size_of::<u32>(),
                        stream,
                    )?;
                }
                // Leases keep these allocations alive until the entire stream
                // submission, including both D2H copies, has completed.
                CompletedOutput::Candidates(CompletedCandidates {
                    keys,
                    ids,
                    layout: selection.layout,
                })
            }
            _ => return Err(DeviceError::Internal("Tile output mode mismatch".into())),
        };
        Ok((output, self.distances))
    }
}

impl IntoFuture for Finish {
    type Output = Result<(CompletedOutput, DeviceBox<[f32]>), DeviceError>;
    type IntoFuture =
        <BoxedDeviceOp<(CompletedOutput, DeviceBox<[f32]>)> as IntoFuture>::IntoFuture;
    fn into_future(self) -> Self::IntoFuture {
        self.boxed().into_future()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cuda_async::device_context::{global_policy, init_device_contexts, with_cuda_context};
    use std::ptr;

    #[test]
    #[ignore = "requires a CUDA GPU"]
    fn cub_returns_chunk_winners_and_fallback_rows() {
        init_device_contexts(0, 1).unwrap();
        let context = with_cuda_context(0, Arc::clone).unwrap();
        let policy = global_policy(0).unwrap();
        let limits = topk::Limits::compiled();
        let nr = 3;
        let mut host = HostOutput::new(&context, 0, Some(nr * (2049 + 8193))).unwrap();
        for (nq, k, expected_stride) in [(8193, 2049, 2049), (8192, 50, 50), (8193, 50, 51)] {
            let selection =
                topk::PreparedSelection::new(nr, nq, k, false, limits, &context).unwrap();
            assert_eq!(selection.layout.ref_stride, expected_stride);
            assert_eq!(selection.layout.query_stride, nr);
            let data: Vec<_> = (0..nr * nq).map(|index| index as f32).collect();
            let stream = policy.next_stream().unwrap();
            let distances = crate::gpu::operations::h2d(&data).sync_on(&stream).unwrap();
            let (output, _) = Finish {
                distances,
                host,
                selection: Some(selection),
            }
            .sync_on(&stream)
            .unwrap();
            let CompletedOutput::Candidates(candidates) = output else {
                panic!("Expected candidates")
            };
            assert_eq!(candidates.layout, selection.layout);
            let CompletedCandidates { keys, ids, .. } = &candidates;
            let mut expected: Vec<_> = (0..k as u32).collect();
            if expected_stride > k {
                expected.push(8192);
            }
            for row in 0..nr {
                let offsets = row * expected_stride..(row + 1) * expected_stride;
                let mut chosen = ids[offsets.clone()].to_vec();
                chosen.sort_unstable();
                assert_eq!(chosen, expected);
                for offset in offsets {
                    assert_eq!(keys[offset], data[row * nq + ids[offset] as usize]);
                }
            }
            let ref_count = nr * expected_stride;
            for column in 0..nq {
                let offsets = ref_count + column * nr..ref_count + (column + 1) * nr;
                let mut chosen = ids[offsets.clone()].to_vec();
                chosen.sort_unstable();
                assert_eq!(chosen, [0, 1, 2]);
                for offset in offsets {
                    assert_eq!(keys[offset], data[ids[offset] as usize * nq + column]);
                }
            }
            host = candidates.into_buffers();
        }
        for (samples, stride) in [(8193, 50), (8194, 51)] {
            let selection =
                topk::PreparedSelection::new(samples, samples, 50, true, limits, &context).unwrap();
            assert_eq!(selection.layout.ref_stride, stride);
            assert_eq!(selection.layout.query_stride, 0);
            let distances: Vec<_> = (0..samples)
                .flat_map(|row| ((row + 1)..samples).map(move |column| (row + column) as f32))
                .collect();
            let stream = policy.next_stream().unwrap();
            let distances = crate::gpu::operations::h2d(&distances)
                .sync_on(&stream)
                .unwrap();
            let host = HostOutput::new(&context, 0, Some(selection.layout.count)).unwrap();
            let (host, _) = Finish {
                distances,
                host,
                selection: Some(selection),
            }
            .sync_on(&stream)
            .unwrap();
            let CompletedOutput::Candidates(CompletedCandidates { keys, ids, layout }) = host
            else {
                panic!("Expected candidates")
            };
            assert_eq!(layout, selection.layout);
            for row in 0..samples {
                let offsets = row * stride..(row + 1) * stride;
                let mut chosen: Vec<_> =
                    ids[offsets.clone()].iter().map(|&id| id as usize).collect();
                chosen.sort_unstable();
                let mut expected: Vec<_> = (0..samples).filter(|&id| id != row).take(50).collect();
                if stride > 50 {
                    expected.push(samples - 1 - usize::from(row == samples - 1));
                }
                assert_eq!(chosen, expected);
                for offset in offsets {
                    assert_eq!(keys[offset], (row + ids[offset] as usize) as f32);
                }
            }
        }
    }

    // Hardware integration coverage for context binding and completion ownership.
    // Selection boundaries are checked separately; this small matrix makes both
    // candidate directions easy to verify while deliberately disturbing contexts.
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
        let limits = topk::Limits::compiled();
        let mut jobs = Vec::new();
        for (id, context) in contexts.iter().enumerate() {
            let foreign = &contexts[(id + 1) % 3];
            let foreign_buffer = PinnedHostBuffer::<u32>::zeroed(foreign, 1).unwrap();
            context.bind_to_thread().unwrap();
            // Foreign resource cleanup changes the current context. Workspace
            // preparation must bind this device again, even for a size query.
            drop(foreign_buffer);
            let selection = topk::PreparedSelection::new(5, 7, 3, false, limits, context).unwrap();
            let mut current = ptr::null_mut();
            unsafe {
                cuda_core::sys::cuCtxGetCurrent(&mut current)
                    .result()
                    .unwrap();
            }
            assert_eq!(current, context.cu_ctx());
            let policy = global_policy(id).unwrap();
            let stream = policy.next_stream().unwrap();
            // Transfer submission must tolerate a different device being current.
            foreign.bind_to_thread().unwrap();
            let distances = crate::gpu::operations::h2d(&data).sync_on(&stream).unwrap();
            let host = HostOutput::new(context, 0, Some(selection.layout.count)).unwrap();
            let job = Finish {
                distances,
                host,
                selection: Some(selection),
            }
            .schedule(policy.as_ref())
            .unwrap();
            // Leave a foreign context current before driving the scheduled work;
            // CUB and raw D2H operations must bind their owning device on execution.
            foreign.bind_to_thread().unwrap();
            jobs.push(job);
        }
        crate::gpu::run_gpu_batches(async {
            let mut pending = tokio::task::JoinSet::new();
            for job in jobs {
                pending.spawn(job);
            }
            while let Some(result) = pending.join_next().await {
                // Completion must retain valid candidate/workspace/transfer storage
                // through both D2H copies before the host reads these results.
                let (host, _distances) = result.unwrap().unwrap();
                let CompletedOutput::Candidates(CompletedCandidates { keys, ids, layout }) = host
                else {
                    panic!("Expected candidates")
                };
                assert_eq!(
                    (
                        layout.k,
                        layout.ref_stride,
                        layout.query_stride,
                        layout.count
                    ),
                    (3, 3, 3, 36)
                );
                // Check reference rows and query columns independently. Sorting
                // IDs avoids imposing an output order on CUB's selected pairs.
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
