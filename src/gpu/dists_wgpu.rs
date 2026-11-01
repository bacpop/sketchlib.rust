use wgpu;
use crate::distances::distance_matrix::{calc_col_idx, calc_query_indices, calc_row_idx, DistType};
use crate::sketch::multisketch::MultiSketch;
use anyhow::{bail, ensure, Context};
use flume::bounded;
use rayon::prelude::*;
use std::collections::{HashMap, VecDeque};
use std::fmt::Write as FmtWrite;
use std::ops::Range;
use std::sync::mpsc::SyncSender;


use wgpu::util::DeviceExt; // This helps you creating and initialising a buffer, either for reading or writing

// Invocations per workgroup that copy the query sketch into workgroup memory: CUDA's
// warp_size (gpu.hpp). Could later be set from the adapter's subgroup size, so that we're adaptive on them. This could enhance performance on the different platforms
// I read that Intel's can go as low as 8, whereas AMD by default (and actually the ONLY valid value for older GPUs, IIRC) uses 64.
const QUERY_LOADERS: u32 = 32;

// CUDA's block size (dist.cu:335-336): used to size the tile ranges, so it always fits,
// and for the interleaving
const DEFAULT_WORKGROUP_SIZE: u32 = 256;
// Smallest workgroup size tile_workgroup_size() may choose, so that the cap on resident
// workgroups per compute unit doesn't limit occupancy (see the notes there)
const MIN_WORKGROUP_SIZE: u32 = 128;
// Whether to choose the workgroup size per tile (tile_workgroup_size()) rather than use
// DEFAULT_WORKGROUP_SIZE for every tile. This could, in principle, be useful but I think other aspects
// like memory usage affect. On my Intel it wasn't worth vs fixed 256 with large indices.
const CHOOSE_WORKGROUP_SIZE_PER_TILE: bool = false;
// Most samples per tile range (the CUDA code's "sketch block"), as its default
// --max-sketch-batch. The device limits can lower it (max_tilerange_size())
const MAX_TILERANGE_SIZE: usize = 10_000;
// Tiles queued on the GPU at once, as the CUDA code's MAX_IN_FLIGHT_TILES. With wgpu's
// single queue and one tile formatted at a time, 2 would already overlap the CPU and the
// GPU; the extra slots absorb tiles of different lengths (diagonal ones have about half the
// pairs), at one largest-tile staging buffer each
const MAX_IN_FLIGHT_TILES: usize = 4;
// Distances per chunk of formatted text, as the CUDA code's DENSE_CHUNK_SIZE
const DENSE_CHUNK_SIZE: usize = 1000;

#[repr(C, align(16))]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct SketchStrides {   // 32 BYTES
    bin_stride: u32,     // 4 bytes
    kmer_stride: u32,    // 4 bytes
    sample_stride: u32,  // 4 bytes
    sketchsize64: u32,   // 4 bytes

    bbits: u32,          // 4 bytes
    _padding: [u32; 3],  // 12 bytes
}

#[repr(C, align(16))]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct DistanceParams {            // 96 BYTES
    ref_strides: SketchStrides,    // 32 bytes
    query_strides: SketchStrides,  // 32 bytes

    ref_n: u32,                   // 4 bytes
    query_n: u32,                 // 4 bytes
    kmer_n: u32,                  // 4 bytes
    dist_n: u32,                  // 4 bytes
    self_comparison: u32,         // 4 bytes
    use_shared: u32,              // 4 bytes

    tilerange_ref_offset: u32,    // 4 bytes: first sample of tilerange_ref
    tilerange_query_offset: u32,  // 4 bytes: first sample of tilerange_query
}

// NOTE on padding: we have 16 bytes WORDS we must align to.
// Type             Size        Alignment
// =======================================
// u32, i32, f32    4 bytes     4 bytes
// vec2<u32>        8 bytes     8 bytes
// vec3<u32>        12 bytes    16 bytes
// vec4<u32>        16 bytes    16 bytes

// Uniform of interleave_sketches.wgsl, for one tile range
#[repr(C)]
#[derive(Clone, Copy, Debug, bytemuck::Pod, bytemuck::Zeroable)]
struct InterleaveParams { // 16 BYTES
    n_samples: u32,       // samples in this tile range
    n_kmers: u32,
    n_bins: u32,          // words per k-mer sketch, all bit planes included
    bin_stride: u32,      // the tile-range size, padded to a multiple of 32
}

// The distance pipeline for one workgroup size, with its bind groups (they belong to the
// pipeline's layout). A diagonal tile reads its one tile range through both sketch bindings
struct DistsPipeline {
    pipeline: wgpu::ComputePipeline,
    diagonal_bg: wgpu::BindGroup,
    off_diagonal_bg: wgpu::BindGroup,
}

/// One pair of tile ranges, computed in one dispatch (the CUDA code's `DistTile`).
/// Rows are tilerange_query and columns tilerange_ref, as in pp-sketchlib's kernel
pub struct DistTile {
    /// Samples of the rows
    pub tilerange_query: Range<usize>,
    /// Samples of the columns
    pub tilerange_ref: Range<usize>,
    /// Both tile ranges are the same one, compared all-vs-all
    pub diagonal: bool,
}

impl DistTile {
    /// Number of distances in this tile
    pub fn n_distances(&self) -> usize {
        let n_ref = self.tilerange_ref.len();
        if self.diagonal {
            n_ref * n_ref.saturating_sub(1) / 2
        } else {
            self.tilerange_query.len() * n_ref
        }
    }

    // Tile-local (query, ref) indices of distance `idx`, as the CUDA code's dist_pair_indices()
    fn dist_pair_indices(&self, idx: usize) -> (usize, usize) {
        let n_ref = self.tilerange_ref.len();
        if self.diagonal {
            let i = calc_row_idx(idx, n_ref);
            (i, calc_col_idx(idx, i, n_ref))
        } else {
            calc_query_indices(idx, n_ref)
        }
    }

    /// Formats this tile's distances in parallel chunks of text, each sent to the writer
    /// as soon as it is ready, in no particular order (the CUDA code's consume(), dense mode)
    pub fn consume(
        &self,
        distances: &[f32],
        names: &[&str],
        text_tx: &SyncSender<String>,
    ) -> anyhow::Result<()> {
        let n_ref = self.tilerange_ref.len();
        distances
            .par_chunks(DENSE_CHUNK_SIZE)
            .enumerate()
            .try_for_each(|(chunk_idx, chunk)| -> anyhow::Result<()> {
                let (mut i, mut j) = self.dist_pair_indices(chunk_idx * DENSE_CHUNK_SIZE);
                let mut text = String::with_capacity(chunk.len() * 24);
                for &distance in chunk {
                    writeln!(
                        text,
                        "{}\t{}\t{distance}",
                        names[self.tilerange_query.start + i],
                        names[self.tilerange_ref.start + j]
                    )?;
                    j += 1;
                    if j == n_ref {
                        i += 1;
                        j = if self.diagonal { i + 1 } else { 0 };
                    }
                }
                text_tx.send(text).context("Distance writer stopped")
            })
    }
}


/// All-vs-all Jaccard distances (or ANI), computed tile by tile as in the CUDA code, with up
/// to MAX_IN_FLIGHT_TILES tiles queued on the GPU. `on_tile` gets every tile's distances as
/// soon as they are read back, in the kernel's dist_idx order within the tile
// TODO: subject to a potential renaming to make everything more uniform and synch with CPU's and CUDA's code, so that it's not a nightmare to support all
pub async fn get_distances_wgpu(
    sketches: &MultiSketch,
    n: usize, // number of samples
    dist_type: &DistType,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
    bin_bits: usize, // BB chosen by the CPU dispatch: BIN_BITS, or LEGACY_BIN_BITS
    mut on_tile: impl FnMut(&DistTile, &[f32]) -> anyhow::Result<()>,
) -> anyhow::Result<()> {
    let &DistType::Jaccard(k_idx, _, ani) = dist_type else {
        bail!("The wgpu backend only supports Jaccard distances for now");
    };
    // No pairs, and wgpu rejects zero-sized buffers
    if n < 2 {
        return Ok(());
    }
    ensure!(
        k_idx < sketches.kmer_lengths().len(),
        "k-mer index {k_idx} out of range"
    );

    // TODO: this could be a parameter? Also, we should expose which GPU to select by index, apart than from PowerPreference, which might be quite useful for browser but a bit less in highly optimised environments.
    let power_pref = wgpu::PowerPreference::None;
    let backends = wgpu::Backends::all();

    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
        backends,
        // This fills the remaining items with the "..THING" and the thing itself.
        ..wgpu::InstanceDescriptor::new_without_display_handle()
    });

    let adapter = instance
        .request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: power_pref, ..Default::default()
        })
        .await
        .context("No GPU adapter found")?;

    // Ask for the hardware limits rather than the WebGPU defaults, so the tile ranges and
    // workgroup sizes reflect what this device can actually do
    let (device, queue) = adapter
        .request_device(&wgpu::DeviceDescriptor {
            label: Some("sketchlib.rust_gpu"),
            required_limits: adapter.limits(),
            ..Default::default()
        })
        .await
        .context("Failed to get GPU device")?;

    log::info!("Using this GPU: {}", adapter.get_info().name);

    // =============================================
    // SETTING UP DATA!
    // (before the pipeline: the overrides depend on the sketch shape)

    // Here we consult the limits of the user's hardware. TODO: use them for QUERY_LOADERS
    // too (the subgroup/warp size)
    let limits = device.limits();
    let max_workgroups = u64::from(limits.max_compute_workgroups_per_dimension);

    // Workgroup sizes to choose from per tile (see tile_workgroup_size()): whole subgroups,
    // from MIN_WORKGROUP_SIZE (or less, on a device whose maximum is lower) up to the
    // device's maximum
    let max_workgroup_size = limits
        .max_compute_workgroup_size_x
        .min(limits.max_compute_invocations_per_workgroup);
    let subgroup_size = adapter.get_info().subgroup_max_size.clamp(1, max_workgroup_size);
    let largest_workgroup_size = max_workgroup_size / subgroup_size * subgroup_size;
    let smallest_workgroup_size = (MIN_WORKGROUP_SIZE.min(max_workgroup_size).div_ceil(subgroup_size)
        * subgroup_size)
        .min(largest_workgroup_size);
    let workgroup_sizes: Vec<u32> = (smallest_workgroup_size..=largest_workgroup_size)
        .step_by(subgroup_size as usize)
        .collect();
    // Always among workgroup_sizes, so every tile has at least one size that fits
    let default_workgroup_size =
        DEFAULT_WORKGROUP_SIZE.min(max_workgroup_size) / subgroup_size * subgroup_size;

    // Sketches as stored, [sample][k-mer][bin]: CUDA's get_strides()
    let (_, n_bins, words_per_sample) = sketches.get_strides();
    let n_kmers = sketches.kmer_lengths().len();
    let all_sketches = sketches.get_sketch_all();

    // Tiles, as the CUDA code's DistTile: tile ranges balanced to at most
    // max_tilerange_size() samples, so that every tile fits this device in one dispatch
    let n_tileranges = n.div_ceil(max_tilerange_size(&limits, default_workgroup_size, words_per_sample)?);
    let tilerange_size = n.div_ceil(n_tileranges);
    let tiles = make_tiles(n, tilerange_size);
    // Per tile, the workgroup size with the fewest idle invocations, or CUDA's for every tile
    let tile_workgroup_sizes = tiles
        .iter()
        .map(|tile| {
            if CHOOSE_WORKGROUP_SIZE_PER_TILE {
                tile_workgroup_size(tile, &workgroup_sizes, max_workgroups)
            } else {
                Ok(default_workgroup_size)
            }
        })
        .collect::<anyhow::Result<Vec<u32>>>()?;
    let mut distinct_sizes = tile_workgroup_sizes.clone();
    distinct_sizes.sort_unstable();
    distinct_sizes.dedup();
    log::info!(
        "{} tile(s), from {n_tileranges} tile range(s) of up to {tilerange_size} samples; workgroup sizes {distinct_sizes:?} (subgroups of up to {subgroup_size}, device maximum {max_workgroup_size})",
        tiles.len()
    );

    // Every tile range is interleaved with the same strides
    let strides = tilerange_strides(sketches, tilerange_size, bin_bits)?;

    // Single k-mer: CUDA would pass `ref + k_idx * kmer_stride`, `kmers + k_idx` and
    // kmer_n = 1. WGSL has no pointer arithmetic, so the offset goes into the binding.
    // We essentially pass as buffer the vector of k values.
    let kmers = [u32::try_from(sketches.kmer_lengths()[k_idx])?];
    let kmer_offset_bytes =
        u64::try_from(k_idx)? * u64::from(strides.kmer_stride) * std::mem::size_of::<u64>() as u64;
    // bin_stride is a multiple of 32 words (256 bytes), which satisfies every adapter
    assert!(
        kmer_offset_bytes % u64::from(limits.min_storage_buffer_offset_alignment) == 0,
        "k-mer offset of {kmer_offset_bytes} bytes is not aligned to {}",
        limits.min_storage_buffer_offset_alignment
    );

    let completeness: Vec<f32> = match completeness_vec {
        Some(cv) => {
            ensure!(
                cv.len() == n,
                "{} completeness values for {n} samples",
                cv.len()
            );
            cv.iter().map(|&c| c as f32).collect()
        }
        // Never read with USE_COMPLETENESS = false, but the binding must exist
        None => vec![1.0],
    };

    // Query sketch in workgroup memory if it fits, as CUDA's check_shared_size()
    // Right now, "query" is also a reference, because we only have all vs all
    let query_words = strides.sketchsize64 * strides.bbits;
    let query_bytes = u64::from(query_words) * std::mem::size_of::<u64>() as u64;
    let use_shared = query_bytes <= u64::from(limits.max_compute_workgroup_storage_size);
    if !use_shared {
        log::warn!(
            "You are using a sketch size that doesn't allow us to use a memory optimisation for the GPGPU code. This might slow down computation on this device. Reduce sketch size to {} or less for better performance",
            64 * u64::from(limits.max_compute_workgroup_storage_size)
                / (u64::from(strides.bbits) * std::mem::size_of::<u64>() as u64)
        );
    }
    log::info!(
        "use_shared={use_shared}: query sketch {query_bytes} bytes, workgroup memory limit {} bytes",
        limits.max_compute_workgroup_storage_size
    );

    // Now we compile the kernel/shader, and prepare to put it as a execution pipeline

    // For embedding at compilation time:
    let shader = device.create_shader_module(wgpu::include_wgsl!("calculate_distances.wgsl"));

    // For runtime recompilation, useful for debugging:
    // let shader_source = std::fs::read_to_string("/home/vrbouza/Projects/Miscellaneous/tests_wgpu/sketchlib.rust/src/gpu/calculate_distances.wgsl").expect("Failed to read WGSL shader");
    // let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
    //         label: Some("calculate_distances"),
    //         source: wgpu::ShaderSource::Wgsl(
    //             shader_source.into()
    //         ),
    //     },
    // );

    // The CUDA code's interleave_sketches kernel
    let interleave_shader =
        device.create_shader_module(wgpu::include_wgsl!("interleave_sketches.wgsl"));
    let interleave_pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("interleave_sketches"),
        layout: None,
        module: &interleave_shader,
        entry_point: None,
        compilation_options: wgpu::PipelineCompilationOptions {
            constants: &[("WORKGROUP_SIZE", f64::from(default_workgroup_size))],
            ..Default::default()
        },
        cache: Default::default(),
    });


    // =============================================
    // We've got the device and other preliminary things set up, we now need to set up the data transfer

    let use_input   = wgpu::BufferUsages::STORAGE  | wgpu::BufferUsages::COPY_DST;
    let use_output  = wgpu::BufferUsages::STORAGE  | wgpu::BufferUsages::COPY_SRC;
    let use_staging = wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST;
    let use_uniform = wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST;
    let create_gpu_buffer = |label: &str, size: u64, usage: wgpu::BufferUsages| {
        device.create_buffer(&wgpu::BufferDescriptor { label: Some(label), size, usage, mapped_at_creation: false })
    };

    let word_bytes = std::mem::size_of::<u64>() as u64;
    // One tile range as stored, and interleaved (samples padded to bin_stride)
    let raw_tilerange_bytes = (tilerange_size * words_per_sample) as u64 * word_bytes;
    let interleaved_tilerange_words = u64::from(strides.kmer_stride) * n_kmers as u64;
    // Sized for the largest tile, and reused by every tile
    let output_bytes = tiles.iter().map(DistTile::n_distances).max().unwrap_or(0) as u64
        * std::mem::size_of::<f32>() as u64;

    // Parameters, written by each tile before its dispatches
    let input_parameters_gpu = create_gpu_buffer("input_parameters_gpu", std::mem::size_of::<DistanceParams>() as u64, use_uniform);
    let interleave_query_params_gpu = create_gpu_buffer("interleave_query_params_gpu", std::mem::size_of::<InterleaveParams>() as u64, use_uniform);
    let interleave_ref_params_gpu = create_gpu_buffer("interleave_ref_params_gpu", std::mem::size_of::<InterleaveParams>() as u64, use_uniform);

    // The tile's two tile ranges as stored (the CUDA code's h2d() of query_bins and
    // ref_bins), and interleaved by interleave_sketches.wgsl
    let tilerange_query_raw_gpu = create_gpu_buffer("tilerange_query_raw_gpu", raw_tilerange_bytes, use_input);
    let tilerange_ref_raw_gpu = create_gpu_buffer("tilerange_ref_raw_gpu", raw_tilerange_bytes, use_input);
    let tilerange_query_gpu = create_gpu_buffer("tilerange_query_gpu", interleaved_tilerange_words * word_bytes, wgpu::BufferUsages::STORAGE);
    let tilerange_ref_gpu = create_gpu_buffer("tilerange_ref_gpu", interleaved_tilerange_words * word_bytes, wgpu::BufferUsages::STORAGE);

    // the k values
    let input_kmers_gpu = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("input_kmers_gpu"),
        contents: bytemuck::cast_slice(&kmers),
        usage: use_input,
    });

    // the completeness, of every sample (as the CUDA code: tiles index it with their offsets)
    let input_completeness_gpu = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("input_completeness_gpu"),
        contents: bytemuck::cast_slice(&completeness),
        usage: use_input,
    });

    let output_dists_gpu = create_gpu_buffer("output_dists_gpu", output_bytes, use_output);
    // One staging buffer per tile in flight, as the CUDA code's pinned host outputs: the CPU
    // reads one while the GPU fills the others
    let n_slots = MAX_IN_FLIGHT_TILES.min(tiles.len());
    let mut free_staging: Vec<wgpu::Buffer> = (0..n_slots)
        .map(|_| create_gpu_buffer("output_dists_staging", output_bytes, use_staging))
        .collect();

    // Interleaving bind groups: raw tile range -> interleaved tile range
    let create_interleave_bind_group = |params: &wgpu::Buffer, raw: &wgpu::Buffer, interleaved: &wgpu::Buffer| {
        device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &interleave_pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: raw.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: interleaved.as_entire_binding() },
            ],
        })
    };
    let interleave_query_bg = create_interleave_bind_group(&interleave_query_params_gpu, &tilerange_query_raw_gpu, &tilerange_query_gpu);
    let interleave_ref_bg = create_interleave_bind_group(&interleave_ref_params_gpu, &tilerange_ref_raw_gpu, &tilerange_ref_gpu);

    // One distance pipeline per workgroup size in use. WGSL fixes the workgroup size when
    // the pipeline is created, as CUDA does at each launch
    let create_dists_pipeline = |workgroup_size: u32| {
        // WebGPU overrides use f64 representation; for bools 0 is false and anything else true
        // This is slightly (I think) better than uniform buffers, that would be the alternative. I've tried to keep equivalence to CUDA's code, unless it could affect seriously performance, to make it easier to mantain.
        let constants = [
            ("WORKGROUP_SIZE", f64::from(workgroup_size)),
            ("SHARED_WORDS", f64::from(if use_shared { query_words } else { 1 })),
            // Every loader must be inside the workgroup, or some words would never be copied
            ("QUERY_LOADERS", f64::from(QUERY_LOADERS.min(workgroup_size))),
            ("ANI", f64::from(u8::from(ani))),
            ("USE_COMPLETENESS", f64::from(u8::from(completeness_vec.is_some()))),
            ("COMPLETENESS_CUTOFF", completeness_cutoff),
        ];
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("calculate_distances"),
            layout: None,
            module: &shader,
            entry_point: None,
            compilation_options: wgpu::PipelineCompilationOptions {
                constants: &constants,
                ..Default::default()
            },
            cache: Default::default(),
        });
        let create_dists_bind_group = |ref_sketches: &wgpu::Buffer| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: None,
                layout: &pipeline.get_bind_group_layout(0),
                entries: &[
                    wgpu::BindGroupEntry { // Parameters
                        binding: 0,
                        resource: input_parameters_gpu.as_entire_binding()
                    },
                    wgpu::BindGroupEntry { // tilerange_ref's sketches
                        binding: 1,
                        // The host-side pointer arithmetic: start at the selected k-mer
                        resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding {
                            buffer: ref_sketches,
                            offset: kmer_offset_bytes, // This is for ONLY Jaccard distances. For when we implement core/acc, we'd need to go over all kmers, but not yet.
                            size: None,
                        }),
                    },
                    wgpu::BindGroupEntry { // output
                        binding: 2,
                        resource: output_dists_gpu.as_entire_binding() },
                    wgpu::BindGroupEntry {// k values
                        binding: 3,
                        resource: input_kmers_gpu.as_entire_binding() },
                    wgpu::BindGroupEntry {// completeness
                        binding: 4,
                        resource: input_completeness_gpu.as_entire_binding() },
                    wgpu::BindGroupEntry { // tilerange_query's sketches, also from the selected k-mer
                        binding: 5,
                        resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding {
                            buffer: &tilerange_query_gpu,
                            offset: kmer_offset_bytes,
                            size: None,
                        }),
                    },
                ],
            })
        };
        DistsPipeline {
            off_diagonal_bg: create_dists_bind_group(&tilerange_ref_gpu),
            diagonal_bg: create_dists_bind_group(&tilerange_query_gpu),
            pipeline,
        }
    };
    let dists_pipelines: HashMap<u32, DistsPipeline> = distinct_sizes
        .iter()
        .map(|&size| (size, create_dists_pipeline(size)))
        .collect();

    // The words of a tile range as stored, the CUDA code's query_bins and ref_bins
    let raw_words = |tilerange: &Range<usize>| {
        &all_sketches[tilerange.start * words_per_sample..tilerange.end * words_per_sample]
    };
    let interleave_params = |tilerange: &Range<usize>| -> anyhow::Result<InterleaveParams> {
        Ok(InterleaveParams {
            n_samples: u32::try_from(tilerange.len())?,
            n_kmers: u32::try_from(n_kmers)?,
            n_bins: u32::try_from(n_bins)?,
            bin_stride: strides.bin_stride,
        })
    };
    // The shader loops over its own invocations, so this stays within the workgroup limit
    let interleave_workgroups = interleaved_tilerange_words
        .div_ceil(u64::from(default_workgroup_size))
        .min(max_workgroups) as u32;

    // The CUDA code's loop over tiles (src/gpu/mod.rs:635-712), without streams: up to
    // n_slots tiles queued on the GPU, each handed to on_tile once read back. The queue runs
    // every write, dispatch and copy in submission order, so the GPU buffers are shared
    let mut tiles_left = tiles.iter().zip(&tile_workgroup_sizes);
    let mut in_flight = VecDeque::new();
    loop {
        // Keep the GPU fed, as the CUDA code's `while !host_outputs.is_empty()`
        while !free_staging.is_empty() {
            let Some((tile, &workgroup_size)) = tiles_left.next() else { break };
            let staging = free_staging.pop().unwrap();

            // Upload the tile ranges (the CUDA code's h2d()); a diagonal tile has just one
            queue.write_buffer(&tilerange_query_raw_gpu, 0, bytemuck::cast_slice(raw_words(&tile.tilerange_query)));
            queue.write_buffer(&interleave_query_params_gpu, 0, bytemuck::bytes_of(&interleave_params(&tile.tilerange_query)?));
            if !tile.diagonal {
                queue.write_buffer(&tilerange_ref_raw_gpu, 0, bytemuck::cast_slice(raw_words(&tile.tilerange_ref)));
                queue.write_buffer(&interleave_ref_params_gpu, 0, bytemuck::bytes_of(&interleave_params(&tile.tilerange_ref)?));
            }

            let n_distances = tile.n_distances();
            let params = DistanceParams {
                ref_strides:     strides,
                query_strides:   strides,
                ref_n:           u32::try_from(tile.tilerange_ref.len())?,
                query_n:         u32::try_from(tile.tilerange_query.len())?,
                kmer_n:          1, // single k-mer, selected by the binding offset
                dist_n:          u32::try_from(n_distances)?,
                self_comparison: u32::from(tile.diagonal), // bools are not supported in structs to parse them to wgsl
                use_shared:      u32::from(use_shared),
                tilerange_ref_offset:   u32::try_from(tile.tilerange_ref.start)?,
                tilerange_query_offset: u32::try_from(tile.tilerange_query.start)?,
            };
            queue.write_buffer(&input_parameters_gpu, 0, bytemuck::bytes_of(&params));

            let mut encoder = device.create_command_encoder(&Default::default());
            {
                let mut pass = encoder.begin_compute_pass(&Default::default());
                // Reorder on the GPU, as the CUDA code's interleave_sketches
                pass.set_pipeline(&interleave_pipeline);
                pass.set_bind_group(0, &interleave_query_bg, &[]);
                pass.dispatch_workgroups(interleave_workgroups, 1, 1);
                if !tile.diagonal {
                    pass.set_bind_group(0, &interleave_ref_bg, &[]);
                    pass.dispatch_workgroups(interleave_workgroups, 1, 1);
                }
                // Distances, with this tile's workgroup size
                let dists = &dists_pipelines[&workgroup_size];
                pass.set_pipeline(&dists.pipeline);
                pass.set_bind_group(0, if tile.diagonal { &dists.diagonal_bg } else { &dists.off_diagonal_bg }, &[]);
                pass.dispatch_workgroups(u32::try_from(tile_workgroups(tile, workgroup_size)?)?, 1, 1);
            } // THis scope is needed because we need to drop pass. We could also drop pass manually
            // Only this tile's distances go back through its staging/intermediate buffer
            let tile_bytes = (n_distances * std::mem::size_of::<f32>()) as u64;
            encoder.copy_buffer_to_buffer(&output_dists_gpu, 0, &staging, 0, tile_bytes);
            let submission = queue.submit([encoder.finish()]);

            // Mapped once the copy is done; checked when this tile's turn comes
            let (tx, rx) = bounded(1);
            staging.map_async(wgpu::MapMode::Read, ..tile_bytes, move |result| {
                // The receiver is gone if an earlier tile failed
                let _ = tx.send(result);
            });
            in_flight.push_back((tile, staging, tile_bytes, submission, rx));
        }

        // The oldest tile: wait for its submission only, so the later ones keep running
        let Some((tile, staging, tile_bytes, submission, rx)) = in_flight.pop_front() else {
            break;
        };
        device.poll(wgpu::PollType::Wait { submission_index: Some(submission), timeout: None })?;
        rx.recv_async().await??;
        {
            let output_data = staging.get_mapped_range(..tile_bytes)?;
            on_tile(tile, bytemuck::cast_slice(&output_data))?;
        }
        // Back to the pool only after consumption, as the CUDA code's host outputs
        staging.unmap();
        free_staging.push(staging);
    }
    Ok(())
}

// Workgroups needed for all-vs-all: each query row of n - q - 1 references gets its
// own workgroups, as in CUDA's getBlockSize()
fn self_dispatch_count(n: u32, workgroup_size: u32) -> u64 {
    (0..n.saturating_sub(1))
        .map(|q| u64::from((n - q - 1).div_ceil(workgroup_size)))
        .sum()
}

// Largest tile range whose tiles each fit this device in one dispatch at workgroup_size. An
// off-diagonal tile is the biggest: size² distances, over ceil(size / workgroup_size)
// workgroups per query row. Each interleaved tile range must also fit a binding and the
// shader's u32 indices
fn max_tilerange_size(
    limits: &wgpu::Limits,
    workgroup_size: u32,
    words_per_sample: usize,
) -> anyhow::Result<usize> {
    let max_binding = (limits.max_storage_buffer_binding_size as u64).min(limits.max_buffer_size);
    let max_workgroups = u64::from(limits.max_compute_workgroups_per_dimension);
    (2..=MAX_TILERANGE_SIZE)
        .rev()
        .find(|&size| {
            let size = size as u64;
            let tilerange_words = size.div_ceil(32) * 32 * words_per_sample as u64;
            size.div_ceil(u64::from(workgroup_size)) * size <= max_workgroups
                && size * size * std::mem::size_of::<f32>() as u64 <= max_binding
                && tilerange_words * std::mem::size_of::<u64>() as u64 <= max_binding
                && tilerange_words <= u64::from(u32::MAX)
        })
        .context("This device can't fit a tile of two samples")
}

// Tile ranges of tilerange_size samples (the last may be smaller). Tiles are the pairs of
// tile ranges i <= j in row-major order, as the CUDA code's tile index; tiles without pairs
// (a tile range of one sample against itself) are dropped, as CUDA does
fn make_tiles(n: usize, tilerange_size: usize) -> Vec<DistTile> {
    let tileranges: Vec<Range<usize>> = (0..n)
        .step_by(tilerange_size)
        .map(|start| start..n.min(start + tilerange_size))
        .collect();
    let mut tiles = Vec::new();
    for (i, tilerange_query) in tileranges.iter().enumerate() {
        for (j, tilerange_ref) in tileranges.iter().enumerate().skip(i) {
            let tile = DistTile {
                tilerange_query: tilerange_query.clone(),
                tilerange_ref: tilerange_ref.clone(),
                diagonal: i == j,
            };
            if tile.n_distances() > 0 {
                tiles.push(tile);
            }
        }
    }
    tiles
}

// Workgroups for one tile: CUDA's getBlockSize() (dist.cu:326-347)
fn tile_workgroups(tile: &DistTile, workgroup_size: u32) -> anyhow::Result<u64> {
    let n_ref = u32::try_from(tile.tilerange_ref.len())?;
    Ok(if tile.diagonal {
        self_dispatch_count(n_ref, workgroup_size)
    } else {
        u64::from(n_ref.div_ceil(workgroup_size)) * tile.tilerange_query.len() as u64
    })
}

// Workgroup size for one tile, among the allowed ones (ascending). Design notes, as neither
// CUDA version does quite this:
//
// - Why per tile: as in pp-sketchlib's kernel, each query row of a tile is covered by whole
//   workgroups, so the end of every row leaves invocations idle: ceil(L / W) * W - L for a
//   row of L references and workgroups of W. pp-sketchlib's getBlockSize() (dist.cu:335-336)
//   only trims W for rows shorter than 256, and the CUDA Rust port always uses 256. Here
//   each tile takes the allowed W with the fewest total invocations, which is the one with
//   the fewest idle invocations along the ref axis.
// - Why idle invocations matter: an idle invocation still reaches every workgroupBarrier(),
//   so it holds its slot on the compute unit until its whole workgroup ends. That slot
//   can't run useful work from another workgroup, so idle invocations lower occupancy.
// - Allowed sizes: multiples of the adapter's largest subgroup size, up to the device's
//   maximum workgroup size. The largest subgroup size is a multiple of whatever size the
//   driver picks, so no subgroup is ever left partly filled (CUDA's warps of 32).
// - Why never below MIN_WORKGROUP_SIZE (128): on its own, "fewest idle invocations" picks
//   tiny workgroups (32 for a diagonal tile of 1,859 samples). But each compute unit also
//   caps how many workgroups it holds at once, whatever their size, and with small
//   workgroups that cap is reached before the unit is full. On NVIDIA (CUDA Programming
//   Guide, compute capabilities), workgroups of 32 reach only 50% occupancy on Hopper and
//   Blackwell (32 blocks x 32 of 2,048 threads per SM) and 33% on Ampere consumer cards
//   (16 x 32 of 1,536). From 128 the cap never binds on any GPU whose limits are documented
//   (Ampere: 16 x 128 = 2,048 >= 1,536). AMD's doesn't bind there either (single-wave
//   workgroups on CDNA; barriers never limit occupancy on RDNA), and Intel and Apple don't
//   publish theirs. wgpu doesn't report these per-unit limits, so the floor is fixed.
// - Ties go to the larger size: fewer workgroups, so fewer copies of the query sketch into
//   workgroup memory, and fewer barriers.
// - Sizes that would need more workgroups than one dispatch allows are skipped. The tile
//   ranges are sized for DEFAULT_WORKGROUP_SIZE (256, CUDA's), so at least one size fits.
// - WGSL fixes the workgroup size when a pipeline is created, so every size in use gets its
//   own pipeline, whereas CUDA sets it at each launch.
// - None of this guarantees full occupancy: registers, and sizes that don't divide the
//   unit's invocation limit (e.g. 3 x 480 = 1,440 of 1,536 on Ada), can still leave slots
//   empty.
// - Only used when CHOOSE_WORKGROUP_SIZE_PER_TILE is on, which it isn't for now (see there).
fn tile_workgroup_size(
    tile: &DistTile,
    workgroup_sizes: &[u32],
    max_workgroups: u64,
) -> anyhow::Result<u32> {
    let mut best: Option<(u64, u32)> = None; // (invocations, workgroup size)
    for &size in workgroup_sizes {
        let workgroups = tile_workgroups(tile, size)?;
        if workgroups > max_workgroups {
            continue;
        }
        let invocations = workgroups * u64::from(size);
        if best.map_or(true, |(fewest, _)| invocations <= fewest) {
            best = Some((invocations, size));
        }
    }
    best.map(|(_, size)| size)
        .context("No workgroup size fits this tile in one dispatch")
}



// The layout of CUDA's flatten_by_samples(), which interleave_sketches.wgsl produces for each tile range on the GPU (as the CUDA code's interleave_sketches kernel).
//
// Originally, the sketches are stored contiguously by sample:
//
//     Sample 0: [A0, A1, A2, A3, ...]
//     Sample 1: [B0, B1, B2, B3, ...]
//     Sample 2: [C0, C1, C2, C3, ...]
//
// The shader transposes them so that samples are the innermost dimension:
//
//     Bin 0: [A0, B0, C0, padding ...]
//     Bin 1: [A1, B1, C1, padding ...]
//     Bin 2: [A2, B2, C2, padding ...]
//     Bin 3: [A3, B3, C3, padding ...]
//
// Each row is padded to a multiple of 32 u64 elements, matching CUDA's
// warpPad() function. The complete layout is [kmer][bin][sample].
//
// WHY?
// In the CUDA distance kernel, each thread in a block compares the same
// query genome against a different reference genome. All threads process
// the same sketch bin at approximately the same time.
//
// By placing the corresponding reference words next to each other in
// memory, neighbouring threads read neighbouring addresses:
//
//     Thread 0 -> A0
//     Thread 1 -> B0
//     Thread 2 -> C0
//
// This allows coalesced memory accesses, reducing the number of memory
// transactions. We preserve this layout in WGSL for compatibility with
// the CUDA implementation, although optimal alignment may differ across
// GPU architectures.
//
// The resulting strides are measured in u64 elements (not bytes):
//
//     sample_stride = 1
//     bin_stride    = the tile-range size, padded to a multiple of 32
//     kmer_stride   = bin_stride * number of words per k-mer sketch
//
// Address of a word:
//     offset = kmer_idx * kmer_stride
//            + bin_idx  * bin_stride
//            + sample_idx * sample_stride

// NOTE TODO: a cool thing would be to make the subgroup/warp size depending on the architecture and make detected, i.e. on the actual GPU you're going to run.

fn tilerange_strides(
    sketches: &MultiSketch,
    tilerange_size: usize,
    bin_bits: usize,
) -> anyhow::Result<SketchStrides> {
    // Sketch size is measured in groups of 64 bins.
    let sketchsize64 = usize::try_from(sketches.sketchsize64)?;
    ensure!(sketchsize64 > 0, "Invalid sketch size");

    // Each k-mer sketch consists of 'sketchsize64 * bbits' u64 words.
    let (_, num_bins, _) = sketches.get_strides();
    ensure!(
        num_bins % sketchsize64 == 0,
        "Sketch word count is not divisible by sketchsize64"
    );
    let bbits = num_bins / sketchsize64;
    ensure!(
        bbits == bin_bits,
        "Sketches have {bbits} bits per bin, but the CPU dispatch expects {bin_bits}"
    );

    // Same stride calculation as CUDA's flatten_by_samples(), for one tile range
    let bin_stride = tilerange_size.div_ceil(32) * 32;
    let kmer_stride = bin_stride
        .checked_mul(num_bins)
        .context("k-mer stride overflow")?;

    // Convert to the u32 fields used by the WGSL uniform.
    Ok(SketchStrides {
        bin_stride: u32::try_from(bin_stride)?,
        kmer_stride: u32::try_from(kmer_stride)?,
        sample_stride: 1,
        sketchsize64: u32::try_from(sketchsize64)?,
        bbits: u32::try_from(bbits)?,
        _padding: [0; 3],
    })
}
