use std::collections::BinaryHeap;
use std::fmt::Write as FmtWrite;
use std::future::Future;
use std::io::Write as IoWrite;
use std::mem;
use std::ops::Range;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, SyncSender};
use std::sync::Arc;

use anyhow::{anyhow, Context, Error};
use rayon::prelude::*;

use cuda_async::device_box::DeviceBox;
use cuda_async::device_context::{init_device_contexts, with_cuda_context};
use cuda_async::device_operation::{self, value, DeviceOperation, Zippable};
use cuda_async::zip;
use cuda_core::device_buffer::DeviceCopy;
use cuda_core::memory::{malloc_async, memcpy_htod_async, memset_d8_async};
use cuda_core::{LaunchConfig1D, PinnedHostBuffer};
use cuda_device::{device, kernel, launch_bounds, launch_contract, thread, DisjointSlice};
use cuda_host::cuda_module;
use tokio::sync::oneshot;
use tokio::task::JoinSet;

use crate::distances::distance_matrix::*;
use crate::distances::jaccard::ani_pois;
use crate::distances::{push_heap, write_jaccard_heap};
use crate::sketch::multisketch::MultiSketch;
use crate::sketch::BIN_BITS;
use crate::utils::get_progress_bar;

mod topk;
use topk::HostOutput;

const BBITS: u64 = BIN_BITS as u64;
const DENSE_CHUNK_SIZE: usize = 1000;

/// Keep the async boundary here; callers do not need a runtime or async API.
fn run_gpu_batches<T>(work: impl Future<Output = Result<T, Error>>) -> Result<T, Error> {
    tokio::runtime::Builder::new_current_thread()
        .build()?
        .block_on(work)
}

/// Decode the same tile-local pair on the device and in the host consumer.
#[device]
fn dist_pair_indices(
    idx: usize,
    n_ref_samples: usize,
    n_query_samples: usize,
    diagonal: bool,
) -> (usize, usize) {
    if diagonal {
        let i = calc_row_idx(idx, n_ref_samples);
        (i, calc_col_idx(idx, i, n_ref_samples))
    } else {
        calc_query_indices(idx, n_query_samples)
    }
}

#[device]
fn calc_jaccard_index(
    sketch1: &[u64],
    sketch2: &[u64],
    sketch1_stride: u64,
    sketch2_stride: u64,
    sketchsize64: u64,
) -> f64 {
    let unionsize = (u64::BITS as u64 * sketchsize64) as f64;
    let mut samebits: u32 = 0;
    for i in 0..sketchsize64 {
        let mut bits: u64 = !0;
        for j in 0..BBITS {
            bits &= !(sketch1[((i * BBITS + j) * sketch1_stride) as usize]
                ^ sketch2[((i * BBITS + j) * sketch2_stride) as usize]);
        }
        samebits += bits.count_ones();
    }
    // Correction for random matches
    let expected_random = unionsize / (1u64 << BBITS) as f64;
    let jaccard_index =
        ((samebits as f64 - expected_random) / (unionsize - expected_random)).clamp(0.0, 1.0);

    jaccard_index
}

#[cuda_module]
mod kernels {
    use super::*;

    #[kernel]
    #[launch_bounds(256)]
    #[launch_contract(domain = 1, block = (256, 1, 1))]
    pub fn dists(
        ref_sketches: &[u64],
        query_sketches: &[u64],
        n_ref_samples: usize,
        n_query_samples: usize,
        ref_strides: (usize, usize, usize), // (bin, k-mer, sample)
        query_strides: (usize, usize, usize), // (bin, k-mer, sample)
        sketchsize64: u64,
        _kmers: &[usize],
        _dist_type: DistType,
        k_vals: Option<(usize, f64)>,
        ani: bool,
        inverse_ani: bool,
        diagonal: bool,
        completeness: &[f32],
        correct_completeness: bool,
        completeness_cutoff: f32,
        sample_offsets: (usize, usize),
        mut dist: DisjointSlice<f32>,
    ) {
        let idx = thread::index_1d();
        let idx_raw = idx.get();

        if let Some(dist_elem) = dist.get_mut(idx) {
            let (i, j) = dist_pair_indices(idx_raw, n_ref_samples, n_query_samples, diagonal);

            if let Some((k_idx, k_f64)) = k_vals {
                let i_idx = i * ref_strides.2 + k_idx * ref_strides.1;
                let j_idx = j * query_strides.2 + k_idx * query_strides.1;

                let mut jaccard_index = calc_jaccard_index(
                    &ref_sketches[i_idx..],
                    &query_sketches[j_idx..],
                    ref_strides.0 as u64,
                    query_strides.0 as u64,
                    sketchsize64,
                );

                if correct_completeness {
                    let c1 = completeness[sample_offsets.0 + i];
                    let c2 = completeness[sample_offsets.1 + j];
                    let product = c1 * c2;
                    // Use the CPU pair-product rule, with f32 correction arithmetic.
                    // Rounding near the cutoff can select a different branch from f64.
                    if product >= completeness_cutoff {
                        jaccard_index = ((jaccard_index as f32 / (product / (c1 + c2 - product)))
                            .min(1.0)) as f64;
                    }
                }

                *dist_elem = if ani {
                    let similarity = ani_pois(jaccard_index, k_f64);
                    if inverse_ani {
                        (1.0_f64 - similarity) as f32
                    } else {
                        similarity as f32
                    }
                } else {
                    (1.0_f64 - jaccard_index) as f32
                };
            } else {
                todo!("path iterating over k-mers and doing regression")
            }
        }
    }

    #[kernel]
    #[launch_bounds(256)]
    #[launch_contract(domain = 1, block = (256, 1, 1))]
    // Host [sample][k-mer][bin] -> device [k-mer][bin][sample].
    // Strides (bin, k-mer, sample): (1, n_bins, n_bins * n_kmers)
    // -> (n_samples, n_samples * n_bins, 1). n_bins includes all bit planes.
    pub fn interleave_sketches(
        sketch_slice: &[u64],
        n_samples: usize,
        n_kmers: usize,
        n_bins: usize,
        mut interleaved: DisjointSlice<u64>,
    ) {
        let idx = thread::index_1d();
        let idx_raw = idx.get();

        // In the output, k-mer outer index, then bin, then sample
        if let Some(interleaved_elem) = interleaved.get_mut(idx) {
            let sample_idx = idx_raw % n_samples;
            let bin_idx = (idx_raw / n_samples) % n_bins;
            let kmer_idx = idx_raw / (n_samples * n_bins);

            let sample_stride_in = n_bins * n_kmers;
            let kmer_stride_in = n_bins;
            *interleaved_elem =
                sketch_slice[sample_idx * sample_stride_in + kmer_idx * kmer_stride_in + bin_idx];
        }
    }
}

/// One pair of sketch blocks, with launch ranges and shared host ownership.
struct DistTile {
    sketches: Arc<PinnedHostBuffer<u64>>,
    ref_samples: Range<usize>,
    query_samples: Range<usize>,
    ref_bins: Range<usize>,
    query_bins: Range<usize>,
    n_kmers: usize,
    n_bins: usize,
    sketchsize64: u64,
    diagonal: bool,
    n_distances: usize,
    interleave_ref_cfg: LaunchConfig1D,
    interleave_query_cfg: LaunchConfig1D,
    distance_cfg: LaunchConfig1D,
}

impl DistTile {
    fn new(
        tile_idx: usize,
        sketches: &MultiSketch,
        max_sketch_batch: usize,
        pinned_sketches: Arc<PinnedHostBuffer<u64>>,
    ) -> Self {
        let n = sketches.number_samples_loaded();
        let n_blocks = n.div_ceil(max_sketch_batch);
        let block_i = calc_row_idx_diag(tile_idx, n_blocks);
        let block_j = calc_col_idx_diag(tile_idx, block_i, n_blocks);
        let start_i = block_i * max_sketch_batch;
        let start_j = block_j * max_sketch_batch;
        let ref_samples = start_i..n.min(start_i + max_sketch_batch);
        let query_samples = start_j..n.min(start_j + max_sketch_batch);
        let n_ref = ref_samples.len();
        let n_query = query_samples.len();
        let (_, n_bins, sample_stride) = sketches.get_strides();
        let ref_bins = ref_samples.start * sample_stride..ref_samples.end * sample_stride;
        let query_bins = query_samples.start * sample_stride..query_samples.end * sample_stride;
        let diagonal = block_i == block_j;
        let n_distances = if diagonal {
            n_ref * (n_ref - 1) / 2
        } else {
            n_ref * n_query
        };

        Self {
            sketches: pinned_sketches,
            interleave_ref_cfg: LaunchConfig1D::new(ref_bins.len().div_ceil(256) as u32, 256, 0),
            interleave_query_cfg: LaunchConfig1D::new(
                query_bins.len().div_ceil(256) as u32,
                256,
                0,
            ),
            distance_cfg: LaunchConfig1D::new(n_distances.div_ceil(256) as u32, 256, 0),
            ref_samples,
            query_samples,
            ref_bins,
            query_bins,
            n_kmers: sketches.kmer_lengths().len(),
            n_bins,
            sketchsize64: sketches.sketchsize64,
            diagonal,
            n_distances,
        }
    }

    /// Build one stream-ordered operation; scheduling and completion stay with
    /// the caller. The returned owners also protect cancelled in-flight work.
    fn pipeline(
        &self,
        module: kernels::LoadedModule,
        kmers: Arc<DeviceBox<[usize]>>,
        dist_type: DistType,
        k_vals: Option<(usize, f64)>,
        ani: bool,
        knn: Option<usize>,
        completeness: Arc<DeviceBox<[f32]>>,
        correct_completeness: bool,
        completeness_cutoff: f32,
        host_output: HostOutput,
    ) -> Result<impl DeviceOperation<Output = (HostOutput, impl Send), IntoFuture: Send> + '_, Error>
    {
        // Prepare contracts before any device work is submitted, so failures
        // propagate normally rather than panicking inside a pipeline closure.
        let prepared_ref = module.prepare_interleave_sketches(self.interleave_ref_cfg)?;
        let prepared_query = module.prepare_interleave_sketches(self.interleave_query_cfg)?;
        let prepared_dist = module.prepare_dists(self.distance_cfg)?;
        let n_ref = self.ref_samples.len();
        let n_query = self.query_samples.len();
        let n_kmers = self.n_kmers;
        let n_bins = self.n_bins;
        let sketchsize64 = self.sketchsize64;
        let diagonal = self.diagonal;
        let n_distances = self.n_distances;
        let sample_offsets = (self.ref_samples.start, self.query_samples.start);
        let query_slice = &self.sketches[self.query_bins.clone()];
        let host_sketches = Arc::clone(&self.sketches);
        let module_owner = module.clone();
        let inverse_ani = knn.is_some();
        let selection = knn
            .map(|k| topk::Shape::new(n_ref, n_query, k, diagonal)?.query())
            .transpose()?;

        // https://nvidia.github.io/cuda-rust/cuda-oxide/async-programming/concurrent-execution.html
        Ok(
            topk::KeepAlive((Arc::clone(&host_sketches), module_owner.clone()))
                .and_then(move |_| {
                    zip!(
                        h2d(&self.sketches[self.ref_bins.clone()]),
                        zeros_u64(self.ref_bins.len())
                    )
                })
                .and_then(move |(input, interleaved)| {
                    module
                        .interleave_sketches_async_owned(
                            &prepared_ref,
                            input,
                            n_ref,
                            n_kmers,
                            n_bins,
                            interleaved,
                        )
                        .and_then(move |(input, interleaved)| value((input, interleaved, module)))
                })
                .and_then(move |(ref_input, ref_interleaved, module)| {
                    zip!(h2d(query_slice), zeros_u64(query_slice.len())).and_then(
                        move |(input, interleaved)| {
                            module
                                .interleave_sketches_async_owned(
                                    &prepared_query,
                                    input,
                                    n_query,
                                    n_kmers,
                                    n_bins,
                                    interleaved,
                                )
                                .and_then(move |(query_input, query_interleaved)| {
                                    value((
                                        ref_input,
                                        ref_interleaved,
                                        query_input,
                                        query_interleaved,
                                        module,
                                    ))
                                })
                        },
                    )
                })
                .and_then(
                    move |(ref_input, ref_interleaved, query_input, query_interleaved, module)| {
                        zeros(n_distances).and_then(move |output| {
                            module
                                .dists_async_owned(
                                    &prepared_dist,
                                    ref_interleaved,
                                    query_interleaved,
                                    n_ref,
                                    n_query,
                                    (n_ref, n_ref * n_bins, 1),
                                    (n_query, n_query * n_bins, 1),
                                    sketchsize64,
                                    kmers,
                                    dist_type,
                                    k_vals,
                                    ani,
                                    inverse_ani,
                                    diagonal,
                                    completeness,
                                    correct_completeness,
                                    completeness_cutoff,
                                    sample_offsets,
                                    output,
                                )
                                .and_then(
                                    move |(
                                        ref_interleaved,
                                        query_interleaved,
                                        kmers,
                                        completeness,
                                        output,
                                    )| {
                                        topk::Finish {
                                            distances: output,
                                            host: host_output,
                                            selection,
                                        }
                                        .and_then(
                                            move |(host, output)| {
                                                // Submission leases independently retain these
                                                // resources through completion or cancellation.
                                                value((
                                                    host,
                                                    (
                                                        ref_input,
                                                        ref_interleaved,
                                                        query_input,
                                                        query_interleaved,
                                                        kmers,
                                                        completeness,
                                                        output,
                                                        host_sketches,
                                                        module_owner,
                                                    ),
                                                ))
                                            },
                                        )
                                    },
                                )
                        })
                    },
                ),
        )
    }

    fn merge_candidates(
        &self,
        host: &HostOutput,
        k: usize,
        heaps: &mut [BinaryHeap<SparseJaccard>],
    ) -> Result<(), Error> {
        let shape = topk::Shape::new(
            self.ref_samples.len(),
            self.query_samples.len(),
            k,
            self.diagonal,
        )?;
        let HostOutput::Candidates { keys, ids } = host else {
            return Err(anyhow!("Expected CUB candidates"));
        };
        let ref_count = self.ref_samples.len() * shape.ref_k;
        for (offset, (&distance, &id)) in keys[..shape.count]
            .iter()
            .zip(&ids[..shape.count])
            .enumerate()
        {
            let (row, neighbour, limit) = if offset < ref_count {
                (
                    self.ref_samples.start + offset / shape.ref_k,
                    self.query_samples.start,
                    self.query_samples.len(),
                )
            } else {
                (
                    self.query_samples.start + (offset - ref_count) / shape.query_k,
                    self.ref_samples.start,
                    self.ref_samples.len(),
                )
            };
            if id as usize >= limit {
                return Err(anyhow!("CUB returned an invalid neighbour ID"));
            }
            push_heap(
                &mut heaps[row],
                SparseJaccard(neighbour + id as usize, distance),
                k,
            );
        }
        Ok(())
    }

    fn consume_dense(
        &self,
        distances: &[f32],
        names: &[&str],
        text_tx: &SyncSender<(String, usize, Option<oneshot::Sender<()>>)>,
        written: oneshot::Sender<()>,
        cancelled: &AtomicBool,
    ) -> Result<(), Error> {
        let n_ref = self.ref_samples.len();
        let n_query = self.query_samples.len();
        distances
            .par_chunks(DENSE_CHUNK_SIZE)
            .enumerate()
            .try_for_each(|(chunk_idx, chunk)| -> Result<(), Error> {
                if cancelled.load(Ordering::Relaxed) {
                    return Err(anyhow!("GPU distance calculation cancelled"));
                }
                let (mut i, mut j) =
                    dist_pair_indices(chunk_idx * DENSE_CHUNK_SIZE, n_ref, n_query, self.diagonal);
                let mut text = String::with_capacity(chunk.len() * 24);
                for &distance in chunk {
                    writeln!(
                        text,
                        "{}\t{}\t{distance}",
                        names[self.ref_samples.start + i],
                        names[self.query_samples.start + j]
                    )?;
                    j += 1;
                    if j == n_query {
                        i += 1;
                        j = if self.diagonal { i + 1 } else { 0 };
                    }
                }
                text_tx
                    .send((text, chunk.len(), None))
                    .context("GPU distance writer stopped")?;
                Ok(())
            })?;
        // All parallel sends have completed, so the marker follows every
        // chunk of this tile even though producers use different threads.
        text_tx
            .send((String::new(), 0, Some(written)))
            .context("GPU distance writer stopped")?;
        Ok(())
    }
}

/// Write GPU self-distances using at most `max_in_flight_tiles` concurrent tiles.
/// Slots remain occupied through CPU consumption; larger limits can increase
/// pinned host and device memory use. This limit does not set the CUDA stream count.
/// `None` streams every unordered pair without a global ordering guarantee;
/// `Some(k)` selects kNN.
/// Completeness correction uses f32, including its pair-product cutoff test;
/// values very close to the cutoff can take a different branch from CPU f64.
/// Direct callers must configure Rayon as for CPU streaming: dense mode needs
/// at least two workers, with one reserved for the writer; kNN needs no reservation.
/// CPU tasks within and across tiles share that pool independently of the tile limit.
/// Sketches must be loaded with `MultiSketch::read_sketch_data_pinned()` or
/// `MultiSketch::read_sketch_data_block_pinned()`; ordinary host storage is rejected.
pub fn self_dists_stream<W: IoWrite + Send>(
    writer: &mut W,
    sketches: &MultiSketch,
    n: usize,
    knn: Option<usize>,
    dist_type: DistType,
    quiet: bool,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
    gpu_device_id: usize,
    max_sketch_batch: usize,
    max_in_flight_tiles: usize,
) -> Result<(), Error> {
    if max_in_flight_tiles == 0 {
        return Err(anyhow!("max_in_flight_tiles must be at least 1"));
    }
    if sketches.is_legacy_format() {
        return Err(anyhow!(
            "GPU dists supports only new-format sketches with BIN_BITS={BIN_BITS}"
        ));
    }

    assert_eq!(n, sketches.number_samples_loaded());
    if n < 256 {
        return Err(anyhow!("GPU self-distances require at least 256 samples"));
    }
    if knn.is_some_and(|knn| knn == 0 || knn >= n) {
        return Err(anyhow!("GPU kNN requires 1 <= knn < n"));
    }
    let (k_vals, ani) = match dist_type {
        DistType::Jaccard(k_idx, k_val, ani) => (Some((k_idx, k_val)), ani),
        DistType::CoreAcc => todo!("GPU core/accessory regression is deferred"),
    };
    let pool_size = rayon::current_num_threads();
    if knn.is_none() && pool_size < 2 {
        return Err(anyhow!("GPU dense streaming requires a Rayon producer and writer; configure the pool with set_threads_with_writer()"));
    }

    let pinned_sketches = Arc::clone(sketches.pinned_sketches().context(
        "GPU distances require pinned sketches; use read_sketch_data_pinned() or read_sketch_data_block_pinned()",
    )?);

    init_device_contexts(gpu_device_id, 1)?;
    // SAFETY: this crate owns the embedded bundle for the kernels above.
    let module = unsafe { kernels::load_async(gpu_device_id)? };
    let kmers_dev = h2d(sketches.kmer_lengths()).arc().sync()?;

    // A disabled correction still passes a valid buffer to the typed launch;
    // its single placeholder value is never indexed by the kernel.
    let completeness: Vec<f32> = completeness_vec.map_or_else(
        || vec![1.0],
        |values| values.iter().map(|&value| value as f32).collect(),
    );
    let completeness_dev = h2d(&completeness).arc().sync()?;
    let correct_completeness = completeness_vec.is_some();
    let mut heaps = if let Some(knn) = knn {
        log::info!("Calculating sparse ref vs ref distances with {knn} nearest neighbours (GPU)");
        // Construct individually: cloning an empty heap loses its reserved capacity.
        // push_heap temporarily needs one spare element before removing the worst.
        (0..n)
            .map(|_| BinaryHeap::with_capacity(knn + 1))
            .collect::<Vec<_>>()
    } else {
        log::info!("Calculating all ref vs ref distances (GPU)");
        Vec::new()
    };
    let names = <DistanceMatrix as Distances>::sketch_names(sketches);
    let n_blocks = n.div_ceil(max_sketch_batch);
    let n_tiles = n_blocks * (n_blocks + 1) / 2;
    let mut tiles = (0..n_tiles)
        .map(|idx| {
            DistTile::new(
                idx,
                sketches,
                max_sketch_batch,
                Arc::clone(&pinned_sketches),
            )
        })
        .filter(|tile| tile.n_distances != 0);
    let initial_tiles: Vec<_> = tiles.by_ref().take(max_in_flight_tiles).collect();
    let n_slots = initial_tiles.len();
    let mut tiles = initial_tiles.into_iter().chain(tiles);
    let block_size = n.min(max_sketch_batch);
    let max_distances =
        (block_size * (block_size - 1) / 2).max(block_size * (n - block_size).min(block_size));
    let max_candidates = knn
        .map(|k| {
            let diagonal = topk::Shape::new(block_size, block_size, k, true)?.count;
            let query_size = (n - block_size).min(block_size);
            let rectangular = if query_size == 0 {
                0
            } else {
                topk::Shape::new(block_size, query_size, k, false)?.count
            };
            Ok::<_, Error>(diagonal.max(rectangular))
        })
        .transpose()?;
    // Allocate before concurrent dispatch: pinned allocations can synchronize.
    let mut host_outputs = with_cuda_context(gpu_device_id, |ctx| {
        (0..n_slots)
            .map(|_| HostOutput::new(ctx, max_distances, max_candidates))
            .collect::<Result<Vec<_>, _>>()
    })??;
    let progress = get_progress_bar(n * (n - 1) / 2, true, quiet);
    let cancelled = AtomicBool::new(false);
    let dense = knn.is_none();
    let consumer_threads = if dense { pool_size - 1 } else { pool_size };

    let mut merge_time = std::time::Duration::ZERO;
    let result = rayon::in_place_scope(|scope| {
        let (text_tx, text_rx) = mpsc::sync_channel::<(
            String, usize, Option<oneshot::Sender<()>>,
        )>(2 * consumer_threads);
        let (writer_done_tx, mut writer_done_rx) = oneshot::channel();
        let (writer_started_tx, writer_started_rx) = oneshot::channel();
        if dense {
            let writer = &mut *writer;
            let progress = &progress;
            let cancelled = &cancelled;
            scope.spawn(move |_| {
                let _ = writer_started_tx.send(());
                let result = catch_unwind(AssertUnwindSafe(|| -> Result<(), Error> {
                    for (text, records, written) in text_rx.iter() {
                        if cancelled.load(Ordering::Relaxed) {
                            return Err(anyhow!("GPU distance calculation cancelled"));
                        }
                        writer.write_all(text.as_bytes())
                            .context("Error writing streamed GPU distance output")?;
                        progress.inc(records as u64);
                        if let Some(written) = written {
                            let _ = written.send(());
                        }
                    }
                    Ok(())
                })).unwrap_or_else(|_| Err(anyhow!("GPU distance writer panicked")));
                // Publish the actual I/O error before receiver drop wakes producers.
                let failed = result.is_err();
                let _ = writer_done_tx.send(result);
                if failed {
                    cancelled.store(true, Ordering::Relaxed);
                }
            });
        } else {
            drop(text_rx);
            drop(writer_done_tx);
            drop(writer_started_tx);
        }

        run_gpu_batches(async {
            if dense {
                // Reserve the writer worker before parallel producers can block
                // on the bounded channel, including with --threads 1.
                writer_started_rx.await.context("GPU distance writer failed to start")?;
            }
            let mut gpu_jobs = JoinSet::<Result<(DistTile, HostOutput), Error>>::new();
            let mut consumers = JoinSet::<Result<HostOutput, Error>>::new();
            let result = async {
                loop {
                    if cancelled.load(Ordering::Relaxed) {
                        (&mut writer_done_rx).await.context("GPU distance writer stopped")??;
                        return Err(anyhow!("GPU distance calculation cancelled"));
                    }
                    // A buffer returns only after consumption, not merely GPU completion.
                    while !host_outputs.is_empty() {
                        let Some(tile) = tiles.next() else { break };
                        let host_output = host_outputs.pop().unwrap();
                        let module = module.clone();
                        let kmers = Arc::clone(&kmers_dev);
                        let completeness = Arc::clone(&completeness_dev);
                        gpu_jobs.spawn(async move {
                            let (host, device_buffers) = tile.pipeline(
                                module, kmers, dist_type, k_vals, ani, knn,
                                completeness, correct_completeness, completeness_cutoff as f32,
                                host_output,
                            )?.await?;
                            // DeviceBox cleanup remains on the GPU dispatch thread.
                            drop(device_buffers);
                            Ok::<_, Error>((tile, host))
                        });
                    }
                    if gpu_jobs.is_empty() && consumers.is_empty() {
                        break;
                    }
                    tokio::select! {
                        biased;
                        writer_result = &mut writer_done_rx, if dense => {
                            writer_result.context("GPU distance writer stopped")??;
                            return Err(anyhow!("GPU distance writer stopped before all tiles were consumed"));
                        }
                        Some(completed) = consumers.join_next(), if !consumers.is_empty() => {
                            let host = completed
                                .context("GPU tile consumer task failed")??;
                            host_outputs.push(host);
                        }
                        Some(completed) = gpu_jobs.join_next(), if !gpu_jobs.is_empty() => {
                            let (tile, host) = completed.context("GPU tile task failed")??;
                            if let Some(k) = knn {
                                let start = std::time::Instant::now();
                                tile.merge_candidates(&host, k, &mut heaps)?;
                                merge_time += start.elapsed();
                                progress.inc(tile.n_distances as u64);
                                host_outputs.push(host);
                                continue;
                            }
                            let (done_tx, done_rx) = oneshot::channel();
                            let (written_tx, written_rx) = oneshot::channel();
                            let text_tx = text_tx.clone();
                            let names = &names;
                            let cancelled = &cancelled;
                            scope.spawn(move |_| {
                                let result = catch_unwind(AssertUnwindSafe(|| {
                                    tile.consume_dense(
                                        &host.dense()[..tile.n_distances], names,
                                        &text_tx, written_tx, cancelled,
                                    )?;
                                    Ok(host)
                                })).unwrap_or_else(|_| Err(anyhow!("GPU tile consumer panicked")));
                                let _ = done_tx.send(result);
                            });
                            consumers.spawn(async move {
                                let result = done_rx.await.context("GPU tile consumer stopped")??;
                                written_rx.await.context("GPU distance writer stopped")?;
                                Ok::<_, Error>(result)
                            });
                        }
                    }
                }
                Ok(())
            }.await;

            if result.is_err() {
                cancelled.store(true, Ordering::Relaxed);
                gpu_jobs.abort_all();
                consumers.abort_all();
                // Join cancellation, not GPU execution: DeviceFuture retains its
                // owned resources behind a completion event until safe reclamation.
                while gpu_jobs.join_next().await.is_some() {}
                while consumers.join_next().await.is_some() {}
            }
            drop(text_tx);
            result?;
            if dense {
                writer_done_rx.await.context("GPU distance writer stopped")??;
            }
            Ok(())
        })
    }).inspect(|_| {
        // Close the distance bar before sparse output moves stderr's cursor.
        progress.finish_and_clear();
    }).inspect_err(|_| {
        progress.abandon();
    }).and_then(|()| {
        if let Some(knn) = knn {
            let mut distances = SparseDistanceMatrix::new(sketches, knn, dist_type);
            if let DistVec::Jaccard(rows) = distances.dists_mut() {
                for (heap, row) in heaps.into_iter().zip(rows.chunks_mut(knn)) {
                    write_jaccard_heap(heap, row, ani);
                }
            }
            log::info!("Writing out in sparse matrix form");
            write!(writer, "{distances}")?;
        }
        Ok(())
    });
    if knn.is_some() {
        log::info!("Host global kNN merging: {:.3}s", merge_time.as_secs_f64());
    }
    result
}

#[allow(dead_code)] // Retain the canonical u64 implementation for reference.
fn h2d_u64(host_data: &[u64]) -> impl DeviceOperation<Output = DeviceBox<[u64]>> + '_ {
    device_operation::with_context(move |ctx| {
        let stream = ctx.get_cuda_stream();
        let n = host_data.len();
        let num_bytes = n * mem::size_of::<u64>();
        unsafe {
            let dptr = malloc_async(stream.cu_stream(), num_bytes).unwrap();
            memcpy_htod_async(dptr, host_data.as_ptr(), num_bytes, stream.cu_stream()).unwrap();
            value(DeviceBox::from_raw_parts(dptr, n, Arc::clone(stream)))
        }
    })
}

/// Generic version of the canonical H2D helper for sketch bins and shared tables.
fn h2d<T: DeviceCopy + Send + Sync + 'static>(
    host_data: &[T],
) -> impl DeviceOperation<Output = DeviceBox<[T]>> + '_ {
    device_operation::with_context(move |ctx| {
        let stream = ctx.get_cuda_stream();
        let n = host_data.len();
        let num_bytes = n * mem::size_of::<T>();
        unsafe {
            let dptr = malloc_async(stream.cu_stream(), num_bytes).unwrap();
            memcpy_htod_async(dptr, host_data.as_ptr(), num_bytes, stream.cu_stream()).unwrap();
            value(DeviceBox::from_raw_parts(dptr, n, Arc::clone(stream)))
        }
    })
}

/// Type-only variant of `zeros()` for interleaved sketch storage.
fn zeros_u64(n: usize) -> impl DeviceOperation<Output = DeviceBox<[u64]>> {
    device_operation::with_context(move |ctx| {
        let stream = ctx.get_cuda_stream();
        let num_bytes = n * mem::size_of::<u64>();
        unsafe {
            let dptr = malloc_async(stream.cu_stream(), num_bytes).unwrap();
            memset_d8_async(dptr, 0, num_bytes, stream.cu_stream()).unwrap();
            value(DeviceBox::from_raw_parts(dptr, n, Arc::clone(stream)))
        }
    })
}

/// Allocates zero-initialized device memory.
/// `ctx` is provided by the scheduler at execution time.
fn zeros(n: usize) -> impl DeviceOperation<Output = DeviceBox<[f32]>> {
    device_operation::with_context(move |ctx| {
        let stream = ctx.get_cuda_stream();
        let num_bytes = n * mem::size_of::<f32>();
        unsafe {
            let dptr = malloc_async(stream.cu_stream(), num_bytes).unwrap();
            memset_d8_async(dptr, 0, num_bytes, stream.cu_stream()).unwrap();
            value(DeviceBox::from_raw_parts(dptr, n, Arc::clone(stream)))
        }
    })
}
