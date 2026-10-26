use std::collections::BinaryHeap;
use std::io::Write as IoWrite;
use std::mem;
use std::ops::Range;
use std::sync::Arc;

use anyhow::{anyhow, Error};

use cuda_async::simt::device_box::DeviceBox;
use cuda_async::simt::device_context::init_device_contexts;
use cuda_async::simt::device_operation::{self, value, DeviceOperation, Zippable};
use cuda_async::zip;
use cuda_core::simt::memory::{
    malloc_async, memcpy_dtoh_async, memcpy_htod_async, memset_d8_async,
};
use cuda_core::{CudaContext, DeviceBuffer, LaunchConfig1D};
use cuda_device::{device, kernel, launch_bounds, launch_contract, thread, DisjointSlice};
use cuda_host::cuda_module;

use crate::distances::distance_matrix::*;
use crate::distances::jaccard::ani_pois;
use crate::distances::{push_heap, write_jaccard_heap};
use crate::sketch::multisketch::MultiSketch;
use crate::sketch::{SketchVec, BIN_BITS};

const BBITS: u64 = BIN_BITS as u64;

// TODO pass async calls through this function
#[cfg(feature = "gpu")]
fn run_gpu_batches() -> Result<(), Error> {
    todo!("asynchronous tile dispatch is deferred")
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
        mut dist: DisjointSlice<f32>,
    ) {
        let idx = thread::index_1d();
        let idx_raw = idx.get();

        if let Some(dist_elem) = dist.get_mut(idx) {
            let (i, j) = dist_pair_indices(idx_raw, n_ref_samples, n_query_samples, diagonal);

            if let Some((k_idx, k_f64)) = k_vals {
                let i_idx = i * ref_strides.2 + k_idx * ref_strides.1;
                let j_idx = j * query_strides.2 + k_idx * query_strides.1;

                let jaccard_index = calc_jaccard_index(
                    &ref_sketches[i_idx..],
                    &query_sketches[j_idx..],
                    ref_strides.0 as u64,
                    query_strides.0 as u64,
                    sketchsize64,
                );

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

/// GPU version of distances::self_dists_all_stream()
pub fn self_dists_all_stream<'a, W: IoWrite + Send>(
    output_file: &mut W,
    sketches: &'a MultiSketch,
    n: usize,
    dist_type: DistType,
    quiet: bool,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
    gpu_device_id: usize,
) -> Result<(), Error> {
    if sketches.is_legacy_format() {
        return Err(anyhow!(
            "GPU dists supports only new-format sketches with BIN_BITS={BIN_BITS}"
        ));
    }

    assert_eq!(n, sketches.number_samples_loaded());
    if n < 256 {
        return Err(anyhow!("GPU kNN requires at least 256 samples"));
    }
    let ctx = CudaContext::new(gpu_device_id)?;
    let stream = ctx.default_stream();

    let mut distances = DistanceMatrix::new(sketches, None, dist_type);
    let k_vals = distances.k_vals();
    let ani = distances.ani();
    let n_distances = distances.n_distances;

    let sketch_dev = DeviceBuffer::from_host(&stream, sketches.get_sketch_all())?;
    let kmers_dev = DeviceBuffer::from_host(&stream, sketches.kmer_lengths())?;

    let mut dist_dev = DeviceBuffer::<f32>::zeroed(&stream, n_distances)?;

    // SAFETY: this package owns the embedded device bundle produced for the
    // kernels module above.
    let module = unsafe { kernels::load(&ctx)? };

    let prepared = module.prepare_dists(LaunchConfig1D::new(
        (n_distances as u32).div_ceil(256),
        256,
        0,
    ))?;

    // TODO interleaving sketches

    module.dists(
        &stream,
        &prepared,
        &sketch_dev,
        &sketch_dev,
        n,
        n,
        sketches.get_strides(),
        sketches.get_strides(),
        sketches.sketchsize64,
        &kmers_dev,
        dist_type,
        k_vals,
        ani,
        false,
        true,
        &mut dist_dev,
    )?;

    let dist_mut = distances.dists_mut();
    *dist_mut = dist_dev.to_host_vec(&stream)?;

    // TODO convert writing to streaming
    log::info!("Writing out in long matrix form");
    write!(output_file, "{distances}")?;

    Ok(())
}

/// One pair of sketch blocks, with launch ranges and shared host ownership.
struct DistTile {
    sketches: Arc<SketchVec>,
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
    fn new(tile_idx: usize, sketches: &MultiSketch, max_sketch_batch: usize) -> Self {
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
            sketches: sketches.get_sketch_all_arc(),
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
    /// the caller. The returned owners keep device allocations alive through sync.
    fn pipeline(
        &self,
        module: kernels::LoadedModule,
        kmers: Arc<DeviceBox<[usize]>>,
        dist_type: DistType,
        k_vals: Option<(usize, f64)>,
        ani: bool,
    ) -> Result<impl DeviceOperation<Output = (Vec<f32>, impl Send)> + '_, Error> {
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
        let query_slice = &self.sketches[self.query_bins.clone()];

        // https://nvidia.github.io/cuda-rust/cuda-oxide/async-programming/concurrent-execution.html
        Ok(zip!(
            h2d(&self.sketches[self.ref_bins.clone()]),
            zeros_u64(self.ref_bins.len())
        )
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
                            true,
                            diagonal,
                            output,
                        )
                        .and_then(move |(ref_interleaved, query_interleaved, kmers, output)| {
                            d2h(output).and_then(move |(host, output)| {
                                // DeviceBox frees on a separate stream. Retain every
                                // allocation in the result until the caller's sync.
                                value((
                                    host,
                                    (
                                        ref_input,
                                        ref_interleaved,
                                        query_input,
                                        query_interleaved,
                                        kmers,
                                        output,
                                    ),
                                ))
                            })
                        })
                })
            },
        ))
    }

    fn consume(&self, distances: Vec<f32>, heaps: &mut [BinaryHeap<SparseJaccard>], knn: usize) {
        // This tile-to-global pair decoding can also feed a future dense-distance
        // streaming writer instead of updating kNN heaps.
        for (idx, distance) in distances.into_iter().enumerate() {
            let (i, j) = dist_pair_indices(
                idx,
                self.ref_samples.len(),
                self.query_samples.len(),
                self.diagonal,
            );
            let i = self.ref_samples.start + i;
            let j = self.query_samples.start + j;
            push_heap(&mut heaps[i], SparseJaccard(j, distance), knn);
            push_heap(&mut heaps[j], SparseJaccard(i, distance), knn);
        }
    }
}

pub fn self_dists_knn<'a>(
    sketches: &'a MultiSketch,
    n: usize,
    knn: usize,
    dist_type: DistType,
    _quiet: bool,
    _completeness_vec: Option<&Vec<f64>>,
    _completeness_cutoff: f64,
    gpu_device_id: usize,
    max_sketch_batch: usize,
) -> Result<SparseDistanceMatrix<'a>, Error> {
    if sketches.is_legacy_format() {
        return Err(anyhow!(
            "GPU dists supports only new-format sketches with BIN_BITS={BIN_BITS}"
        ));
    }

    assert_eq!(n, sketches.number_samples_loaded());
    if n < 256 {
        return Err(anyhow!("GPU kNN requires at least 256 samples"));
    }
    if knn >= n {
        return Err(anyhow!("GPU kNN requires knn to be smaller than n"));
    }
    if matches!(dist_type, DistType::CoreAcc) {
        todo!("GPU core/accessory regression is deferred");
    }

    init_device_contexts(gpu_device_id, 1)?;
    // SAFETY: this crate owns the embedded bundle for the kernels above.
    let module = unsafe { kernels::load_async(gpu_device_id)? };
    let kmers_dev = h2d_usize(sketches.kmer_lengths()).arc().sync()?;

    let mut distances = SparseDistanceMatrix::new(sketches, knn, dist_type);
    let k_vals = distances.k_vals();
    let ani = distances.ani();
    let mut heaps = vec![BinaryHeap::with_capacity(knn + 1); n];
    let n_blocks = n.div_ceil(max_sketch_batch);
    let n_tiles = n_blocks * (n_blocks + 1) / 2;

    for tile_idx in 0..n_tiles {
        let tile = DistTile::new(tile_idx, sketches, max_sketch_batch);
        if tile.n_distances == 0 {
            continue;
        }
        let (result, _device_buffers) = tile
            .pipeline(
                module.clone(),
                Arc::clone(&kmers_dev),
                dist_type,
                k_vals,
                ani,
            )?
            .sync()?;
        tile.consume(result, &mut heaps, knn);
    }

    if let DistVec::Jaccard(rows) = distances.dists_mut() {
        for (heap, row) in heaps.into_iter().zip(rows.chunks_mut(knn)) {
            write_jaccard_heap(heap, row, ani);
        }
    }
    Ok(distances)
}

fn h2d(host_data: &[u64]) -> impl DeviceOperation<Output = DeviceBox<[u64]>> + '_ {
    device_operation::with_context(move |ctx| {
        let stream = ctx.get_cuda_stream();
        let n = host_data.len();
        let num_bytes = n * mem::size_of::<u64>();
        unsafe {
            let dptr = malloc_async(stream.cu_stream(), num_bytes).unwrap();
            memcpy_htod_async(dptr, host_data.as_ptr(), num_bytes, stream.cu_stream()).unwrap();
            value(DeviceBox::from_raw_parts(dptr, n, ctx.get_device_id()))
        }
    })
}

fn h2d_usize(host_data: &[usize]) -> impl DeviceOperation<Output = DeviceBox<[usize]>> + '_ {
    device_operation::with_context(move |ctx| {
        let stream = ctx.get_cuda_stream();
        let n = host_data.len();
        let num_bytes = n * mem::size_of::<usize>();
        unsafe {
            let dptr = malloc_async(stream.cu_stream(), num_bytes).unwrap();
            memcpy_htod_async(dptr, host_data.as_ptr(), num_bytes, stream.cu_stream()).unwrap();
            value(DeviceBox::from_raw_parts(dptr, n, ctx.get_device_id()))
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
            value(DeviceBox::from_raw_parts(dptr, n, ctx.get_device_id()))
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
            value(DeviceBox::from_raw_parts(dptr, n, ctx.get_device_id()))
        }
    })
}

/// Copies device memory back to host (D2H).
/// `ctx` is provided by the scheduler at execution time.
fn d2h(dev: DeviceBox<[f32]>) -> impl DeviceOperation<Output = (Vec<f32>, DeviceBox<[f32]>)> {
    device_operation::with_context(move |ctx| {
        let stream = ctx.get_cuda_stream();
        let n = dev.len();
        let num_bytes = n * mem::size_of::<f32>();
        let mut host = vec![0.0f32; n];
        unsafe {
            memcpy_dtoh_async(
                host.as_mut_ptr(),
                dev.cu_deviceptr(),
                num_bytes,
                stream.cu_stream(),
            )
            .unwrap();
        }
        // Keep dev in the output until the caller has synchronised the stream.
        value((host, dev))
    })
}
