use std::io::Write as IoWrite;

use anyhow::Error;

use cuda_core::{CudaContext, DeviceBuffer, LaunchConfig1D};
use cuda_device::{DisjointSlice, kernel, launch_bounds, launch_contract, thread, device};
use cuda_host::cuda_module;
use cuda_async::simt::device_context::init_device_contexts;

use crate::sketch::multisketch::MultiSketch;
use crate::sketch::BIN_BITS;
use crate::distances::distance_matrix::*;
use crate::distances::jaccard::ani_pois;

const BBITS: u64 = BIN_BITS as u64;

// TODO pass async calls through this function
#[cfg(feature = "gpu")]
fn run_gpu_batches() -> Result<(), Error> {
    let runtime = tokio::runtime::Builder::new_current_thread()
        .enable_all()
        .build()?;

    runtime.block_on(async {
        let batches = prepare_batches();

        let results = process_batches(batches).await?;

        consume_results(results);

        Ok(())
    })
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
            bits &= !(sketch1[((i * BBITS + j) * sketch1_stride) as usize] ^ sketch2[((i * BBITS + j) * sketch2_stride) as usize]);
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
        sketches: &[u64],
        n_ref_samples: usize,
        n_query_samples: usize,
        ref_strides: (usize, usize, usize), // (bin, k-mer, sample)
        query_strides: (usize, usize, usize), // (bin, k-mer, sample)
        sketchsize64: u64,
        _kmers: &[usize],
        _dist_type: DistType,
        k_vals: Option<(usize, f64)>,
        ani: bool,
        diagonal: bool,
        mut dist: DisjointSlice<f32>
    ) {
        let idx = thread::index_1d();
        let idx_raw = idx.get();

        if let Some(dist_elem) = dist.get_mut(idx) {
            let (i, j) = if diagonal {
                (calc_row_idx(idx_raw, n_samples), calc_col_idx(idx_raw, i, n_samples))
            } else {
                calc_query_indices(idx_raw, n_query_samples)
            };

            if let Some((k_idx, k_f64)) = k_vals {
                let i_idx = i * ref_strides.2 + k_idx * ref_strides.1;
                let j_idx = j * query_strides.2 + k_idx * query_strides.1;

                let jaccard_index = calc_jaccard_index(
                                            &sketches[i_idx..],
                                            &sketches[j_idx..],
                                            ref_strides.0 as u64,
                                            query_strides.0 as u64,
                                            sketchsize64);

                *dist_elem = if ani {
                    ani_pois(jaccard_index, k_f64) as f32
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
    // Interleaves 3D array bin, kmer, sample to sample, bin, k-mer (in order of fastest moving index)
    // i.e. [1, n_bins, n_bins * n_kmers] to [n_samples, n_samples * n_bins , 1]
    // TODO check this logic
    pub fn interleave_sketches(
        sketch_slice: &[u64],
        n_samples: usize,
        n_kmers: usize,
        n_bins: usize,
        mut interleaved: DisjointSlice<u64>
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
            *interleaved_elem = sketch_slice
                [sample_idx * sample_stride_in + kmer_idx * kmer_stride_in + bin_idx];
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
    if knn >= n {
        return Err(anyhow!("GPU kNN requires knn to be smaller than n"));
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

    let prepared = module.prepare_dists(
        LaunchConfig1D::new((n_distances as u32).div_ceil(256), 256, 0)
    )?;

    // TODO interleaving sketches

    module.dists(
        &stream,
        &prepared,
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

pub fn self_dists_knn<'a>(
    sketches: &'a MultiSketch,
    n: usize,
    knn: usize,
    dist_type: DistType,
    quiet: bool,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
    gpu_device_id: usize,
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

    init_device_contexts(gpu_device_id, 1)?;
    let module = unsafe { kernels::load_async(0)? };

    let mut sp_distances = SparseDistanceMatrix::new(sketches, knn, dist_type);
    let k_vals = sp_distances.k_vals();
    let ani = sp_distances.ani();

    let (_, _, n_sketch_bins) = sketches.get_strides();

    let kmers_dev = DeviceBuffer::from_host(&stream, sketches.kmer_lengths())?;

    ///////////////

    // TODO 100k sketches per batch
    let num_linear_blocks = 4;
    let num_batches = num_linear_blocks * (num_linear_blocks + 1) / 2;

    let samples_per_block = sketches.number_samples_loaded() / num_linear_blocks;
    let final_block_size = sketches.number_samples_loaded() % num_linear_blocks;

    let mut handles = vec![];

    for batch_idx in 0..num_batches {
        let module = module.clone();

        // TODO doesn't currently work for diagonal blocks
        let block_i = calc_row_idx(batch_idx, num_linear_blocks);
        let block_j = calc_col_idx(batch_idx, block_i, num_linear_blocks);

        let start_i = block_i * samples_per_block;
        let start_j = block_j * samples_per_block;

        let end_i = if block_i < num_linear_blocks - 1 {
            start_i + samples_per_block - 1
        } else {
            start_i + final_block_size - 1
        };

        let end_j = if block_j < num_linear_blocks - 1 {
            start_j + samples_per_block - 1
        } else {
            start_j + final_block_size - 1
        };

        let n_dists_in_block = if block_i == num_linear_blocks - 1 && block_j == num_linear_blocks - 1 {
            final_block_size * (final_block_size - 1) / 2
        } else if block_i == block_j {
            samples_per_block * (samples_per_block - 1) / 2
        //} else if block_j == num_linear_blocks - 1 {
            //0
        } else {
            (end_i - start_i + 1) * (end_j - start_j + 1)
        };

        let ref_sketch_slice = sketches.get_sketch_all()[start_i * n_sketch_bins..end_i * n_sketch_bins];
        let query_sketch_slice = sketches.get_sketch_all()[start_j * n_sketch_bins..end_j * n_sketch_bins];

        let interleave_query_cfg = LaunchConfig::for_num_elems(ref_sketch_slice.len());
        let interleave_ref_cfg = LaunchConfig::for_num_elems(query_sketch_slice.len());
        let distance_cfg = LaunchConfig::for_num_elems(n_dists_in_block);

        // Example code for the pipeline:
        // https://nvidia.github.io/cuda-rust/cuda-oxide/async-programming/concurrent-execution.html
        // TODO: do we need PinnedHostBuffer<T> for the sketch strides we will copy?

        // TODO pipeline steps:
        let pipeline =
            // 1) H2D ref sketch
            zip!(h2d(&ref_sketch_slice), zeros(ref_sketch_slice.len()))
            // 2) Interleave ref sketches
            .and_then(
                move |(input, interleaved): (
                    DeviceBox<[u64]>,
                    DeviceBox<[u64]>,
                )| {
                    let launch = unsafe {
                        module.interleave_sketches_async_owned(
                            input,
                            end_i - start_i + 1,
                            sketches.kmer_lengths(),
                            1,
                            interleaved,
                        )
                    }
                    .expect("Failed to build interleaved launch");

                    // TODO: unsure how to pass interleaved onto next fn
                    launch.and_then(move |(_input, _w0, hidden)| {
                        value((hidden.into_buffer(), output, w1, module))
                    })
                },
            )
            // 3) H2D query sketch
            .and_then(zip!(h2d(&query_sketch_slice), zeros(query_sketch_slice.len())))
            // 4) Interleave query sketches
            .and_then(
                move |(input, interleaved): (
                    DeviceBox<[u64]>,
                    DeviceBox<[u64]>,
                )| {
                    let launch = unsafe {
                        module.interleave_sketches_async_owned(
                            input,
                            end_j - start_j + 1,
                            sketches.kmer_lengths(),
                            1,
                            interleaved,
                        )
                    }
                    .expect("Failed to build interleaved launch");

                    // TODO: unsure how to pass interleaved onto next fn
                    launch.and_then(move |(_input, _w0, hidden)| {
                        value((hidden.into_buffer(), output, w1, module))
                    })
                },
            )
            // 5) Calculate distance kernel
            .and_then(
                // TODO call dists, with correct strides, and correct
            )
            // 6) Copy distances back to the host
            .and_then(d2h);

        // Spawn the pipeline as a Tokio task. The DeviceOperation is
        // converted to a Future via IntoFuture and polled by the runtime.
        handles.push(tokio::spawn(pipeline.into_future()));
    }



    // Await all batches and verify results.
    let knn_heaps = vec![BinaryHeap::with_capacity(knn + 1); sketches.number_samples_loaded()];
    log::info!(
        "Launched {} batches concurrently, awaiting results...\n",
        num_batches
    );
    for (block_idx, handle) in handles.into_iter().enumerate() {
        let result: Vec<f32> = handle.await.expect("Tokio task panicked")?;
        // TODO In consumer push to heaps for (i, j) and (j, i)
        // - calculate the i and j, given block_idx, and index within the block
        // - call push_heap() twice, for heap[i] idx j, and heap[j] idx i
    }
    // TODO call l398-407 from distances/mod.rs (self_dists_knn_generic()) to move
    // heaps into sparse matrix output
    // (move this code into a function called here, and in the original place)
    todo!()
}

fn h2d(host_data: &[u64]) -> impl DeviceOperation<Output = DeviceBox<[u64]>> {
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
fn d2h(dev: DeviceBox<[f32]>) -> impl DeviceOperation<Output = Vec<f32>> {
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
        // dev is kept alive in this closure until the stream synchronizes
        let _ = &dev;
        value(host)
    })
}
