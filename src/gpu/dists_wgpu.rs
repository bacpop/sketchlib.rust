use wgpu;
use crate::distances::distance_matrix::DistType;
use crate::sketch::multisketch::MultiSketch;
use anyhow::{bail, ensure, Context};
use flume::bounded;


use wgpu::util::DeviceExt; // This helps you creating and initialising a buffer, either for reading or writing

// Invocations per workgroup that copy the query sketch into workgroup memory: CUDA's
// warp_size (gpu.hpp). Could later be set from the adapter's subgroup size, so that we're adaptive on them. This could enhance performance on the different platforms
// I read that Intel's can go as low as 8, whereas AMD by default (and actually the ONLY valid value for older GPUs, IIRC) uses 64.
const QUERY_LOADERS: u32 = 32;

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

    _padding: [u32; 2],           // 8 bytes
}

// NOTE on padding: we have 16 bytes WORDS we must align to.
// Type             Size        Alignment
// =======================================
// u32, i32, f32    4 bytes     4 bytes
// vec2<u32>        8 bytes     8 bytes
// vec3<u32>        12 bytes    16 bytes
// vec4<u32>        16 bytes    16 bytes


/// All-vs-all Jaccard distances (or ANI), in condensed upper-triangle order (the same as `square_to_condensed`)
// TODO: subject to a potential renaming to make everything more uniform and synch with CPU's and CUDA's code, so that it's not a nightmare to support all
pub async fn get_distances_wgpu(
    sketches: &MultiSketch,
    n: usize, // number of samples
    dist_type: &DistType,
    completeness_vec: Option<&Vec<f64>>,
    completeness_cutoff: f64,
    bin_bits: usize, // BB chosen by the CPU dispatch: BIN_BITS, or LEGACY_BIN_BITS
) -> anyhow::Result<Vec<f32>> {
    let &DistType::Jaccard(k_idx, _, ani) = dist_type else {
        bail!("The wgpu backend only supports Jaccard distances for now");
    };
    let n = u32::try_from(n)?;
    // No pairs, and wgpu rejects zero-sized buffers
    if n < 2 {
        return Ok(Vec::new());
    }
    let n_pairs = u64::from(n) * u64::from(n - 1) / 2;
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

    // Ask for the hardware limits rather than the WebGPU defaults, so check_limits()
    // reflects what this device can actually do
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

    let (flat_sketches, strides) = flatten_by_samples(sketches)?;
    ensure!(
        strides.bbits as usize == bin_bits,
        "Sketches have {} bits per bin, but the CPU dispatch expects {bin_bits}",
        strides.bbits
    );

    // Single k-mer: CUDA would pass `ref + k_idx * kmer_stride`, `kmers + k_idx` and
    // kmer_n = 1. WGSL has no pointer arithmetic, so the offset goes into the binding.
    // We essentially pass as buffer the vector of k values.
    let kmers = [u32::try_from(sketches.kmer_lengths()[k_idx])?];
    let kmer_offset_bytes =
        u64::try_from(k_idx)? * u64::from(strides.kmer_stride) * std::mem::size_of::<u64>() as u64;

    let completeness: Vec<f32> = match completeness_vec {
        Some(cv) => {
            ensure!(
                cv.len() == n as usize,
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
    // Here we consult the limits of the user's hardware. TODO: use this to update dynamically the
    // subgroup/warp size value, as well as to optimise the dispatches.
    let limits = device.limits();
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

    // Main parameters struct, passed via uniform buffer
    let params = DistanceParams {
        ref_strides:    strides,
        query_strides:  strides,
        ref_n:          n,
        query_n:        n,
        kmer_n:         1, // single k-mer, selected by the binding offset
        dist_n:         u32::try_from(n_pairs)?,
        self_comparison: 1, // bools are not supported in structs to parse them to wgsl
        use_shared:     u32::from(use_shared),
        _padding:        [0; 2],
    };

    // We specified 256 threads per workgroup in the shader, so we need to compute how many
    // TODO: same workgroup/block size adaptation that is done in the CUDA version, AND make it adaptive depending on device limits
    let workgroup_size = 256u32;

    // This is for allvsall
    let num_dispatches = self_dispatch_count(n, workgroup_size);
    let input_bytes = std::mem::size_of_val(flat_sketches.as_slice()) as u64;
    let output_bytes = n_pairs * std::mem::size_of::<f32>() as u64;
    check_limits(&limits, n, num_dispatches, input_bytes, output_bytes, flat_sketches.len());
    // bin_stride is a multiple of 32 words (256 bytes), which satisfies every adapter
    assert!(
        kmer_offset_bytes % u64::from(limits.min_storage_buffer_offset_alignment) == 0,
        "k-mer offset of {kmer_offset_bytes} bytes is not aligned to {}",
        limits.min_storage_buffer_offset_alignment
    );

    // OLD: for queryvsall
    // workgroups we need to dispatch. This code is equivalent to getBlockSize() from the CUDA version. This is for query vs ref, not what we're doing now
    // let ref_samples = sketches.number_samples_loaded();
    // let workgroups_per_query = (ref_samples + workgroup_size - 1) / workgroup_size;
    // let num_dispatches = (workgroups_per_query * ref_samples) as u32;

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

    // WebGPU overrides use f64 representation; for bools 0 is false and anything else true
    // This is slightly (I think) better than uniform buffers, that would be the alternative. I've tried to keep equivalence to CUDA's code, unless it could affect seriously performance, to make it easier to mantain.
    let constants = [
        ("WORKGROUP_SIZE", f64::from(workgroup_size)),
        ("SHARED_WORDS", f64::from(if use_shared { query_words } else { 1 })),
        ("QUERY_LOADERS", f64::from(QUERY_LOADERS)),
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


    // =============================================
    // We've got the device and other preliminary things set up, we now need to set up the data transfer

    let use_input   = wgpu::BufferUsages::STORAGE  | wgpu::BufferUsages::COPY_DST;
    let use_output  = wgpu::BufferUsages::STORAGE  | wgpu::BufferUsages::COPY_SRC;
    let use_staging = wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST;
    let use_uniform = wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST;

    // parameters
    let input_parameters_gpu = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("input_parameters_gpu"),
        contents: bytemuck::bytes_of(&params),
        usage: use_uniform,
    });

    // the sketches
    let input_sketches_gpu = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("input_sketches_gpu"),
        // NOTE: big-endian systems will INVERT the order of how this info is passed to the WGSL kernel, but this happens in weird places, so let's leave it for the time being as it is...
        contents: bytemuck::cast_slice(&flat_sketches),
        usage: use_input,
    });

    // the k values
    let input_kmers_gpu = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("input_kmers_gpu"),
        contents: bytemuck::cast_slice(&kmers),
        usage: use_input,
    });

    // the completeness
    let input_completeness_gpu = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("input_completeness_gpu"),
        contents: bytemuck::cast_slice(&completeness),
        usage: use_input,
    });

    // For outputs, don't use the init one
    let output_dists_gpu = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("output_dists_gpu"),
        size: output_bytes,
        usage: use_output,
        mapped_at_creation: false,
    });

    // We need an intermediate buffer for "mapping" (bringing back to memory) the data to CPU
    let output_dists_staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("output_dists_staging"),
        size: output_bytes,
        usage: use_staging,
        mapped_at_creation: false,
    });

    // println!("{:?}", sketches.get_sketch_slice(0, 0)); // TEST

    // Now, we create the bindgroup: this links one shader with the buffers it needs
    let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label:   None,
        layout:  &pipeline.get_bind_group_layout(0),
        entries: &[
            wgpu::BindGroupEntry { // Parameters
                binding: 0,
                resource: input_parameters_gpu.as_entire_binding()
            },
            wgpu::BindGroupEntry { // input (sketches)
                binding: 1,
                // The host-side pointer arithmetic: start at the selected k-mer
                resource: wgpu::BindingResource::Buffer(wgpu::BufferBinding {
                    buffer: &input_sketches_gpu,
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
        ],
    });


    let mut encoder = device.create_command_encoder(&Default::default());
    {
        let mut pass = encoder.begin_compute_pass(&Default::default());
        pass.set_pipeline(&pipeline);
        pass.set_bind_group(0, &bg, &[]);
        // check_limits() guarantees this fits in u32
        pass.dispatch_workgroups(num_dispatches as u32, 1, 1);
    } // THis scope is needed because we need to drop pass. We could also drop pass manually
    // This is important and allows to get back the results via the staging/intermediate buffer
    encoder.copy_buffer_to_buffer(&output_dists_gpu, 0, &output_dists_staging, 0, output_dists_gpu.size());

    queue.submit([encoder.finish()]);
    let mut output_data : Vec<f32> = Vec::new();

    { // And this scope is also important because of the mapping/unmapping
        // Now, the actual execution
        // The mapping process is async, so we'll need to create a channel to get
        // the success flag for our mapping
        let (tx, rx) = bounded(1);

        // We send the success or failure of our mapping via a callback
        output_dists_staging.map_async(wgpu::MapMode::Read, .., move |result| {
            tx.send(result).unwrap()
        });

        // The callback we submitted to map async will only get called after the
        // device is polled or the queue submitted
        device.poll(wgpu::PollType::wait_indefinitely())?;

        // We check if the mapping was successful here, and we wait for the end of the running
        rx.recv_async().await??;

        // We then get the bytes that were stored in the buffer
        let output_data_ref = &output_dists_staging.get_mapped_range(..)?;
        output_data.extend_from_slice(bytemuck::cast_slice(output_data_ref));

        // Now we have the data on the CPU we can do what ever we want to with it
        // assert_eq!(&input_data, bytemuck::cast_slice(&output_data));
    }

    output_dists_staging.unmap();

    // println!("{:?}", output_data);

    debug_assert_eq!(output_data.len() as u64, n_pairs);
    Ok(output_data)
}

// Workgroups needed for all-vs-all: each query row of n - q - 1 references gets its
// own workgroups, as in CUDA's getBlockSize()
fn self_dispatch_count(n: u32, workgroup_size: u32) -> u64 {
    (0..n.saturating_sub(1))
        .map(|q| u64::from((n - q - 1).div_ceil(workgroup_size)))
        .sum()
}

// Host-side guards for limits that the GPU would otherwise hit with a validation error,
// or silently (out-of-range accesses dropped or clamped, u32 wrap-around)
fn check_limits(
    limits: &wgpu::Limits,
    n: u32,
    num_dispatches: u64,
    input_bytes: u64,
    output_bytes: u64,
    input_words: usize,
) {
    assert!(
        num_dispatches <= u64::from(limits.max_compute_workgroups_per_dimension),
        "{n} samples need {num_dispatches} workgroups, above this device's limit of {}",
        limits.max_compute_workgroups_per_dimension
    );

    let max_binding = (limits.max_storage_buffer_binding_size as u64).min(limits.max_buffer_size);
    assert!(
        input_bytes <= max_binding,
        "Flattened sketches take {input_bytes} bytes, above this device's buffer limit of {max_binding}"
    );
    assert!(
        output_bytes <= max_binding,
        "Distances take {output_bytes} bytes, above this device's buffer limit of {max_binding}"
    );

    // The shader computes q * (2n - q - 1) in u32; it peaks at the last query, q = n - 2
    assert!(
        u64::from(n - 2) * (u64::from(n) + 1) <= u64::from(u32::MAX),
        "{n} samples overflow the shader's u32 condensed-index arithmetic"
    );
    // Every sketch word is addressed with a u32 index
    assert!(
        u32::try_from(input_words).is_ok(),
        "Flattened sketches have {input_words} words, beyond the shader's u32 indexing"
    );
}



// Reorder the sketches to reproduce CUDA's flatten_by_samples() layout.
//
// Originally, the sketches are stored contiguously by sample:
//
//     Sample 0: [A0, A1, A2, A3, ...]
//     Sample 1: [B0, B1, B2, B3, ...]
//     Sample 2: [C0, C1, C2, C3, ...]
//
// We transpose them so that samples are the innermost dimension:
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
//     bin_stride    = number of samples padded to a multiple of 32
//     kmer_stride   = bin_stride * number of words per k-mer sketch
//
// Address of a word:
//     offset = kmer_idx * kmer_stride
//            + bin_idx  * bin_stride
//            + sample_idx * sample_stride

// NOTE TODO: a cool thing would be to make the subgroup/warp size depending on the architecture and make detected, i.e. on the actual GPU you're going to run.

fn flatten_by_samples(
    sketches: &MultiSketch,
) -> anyhow::Result<(Vec<u64>, SketchStrides)> {
    use anyhow::{ensure, Context};

    let n_samples = sketches.number_samples_loaded();
    let n_kmers = sketches.kmer_lengths().len();

    ensure!(n_samples > 0, "No sketches loaded");
    ensure!(n_kmers > 0, "No k-mer lengths available");

    // Sketch size is measured in groups of 64 bins.
    let sketchsize64 = usize::try_from(sketches.sketchsize64)?;

    // Each k-mer sketch consists of 'sketchsize64 * bbits' u64 words.
    let num_bins = sketches.get_sketch_slice(0, 0).len();

    ensure!(sketchsize64 > 0, "Invalid sketch size");
    ensure!(
        num_bins % sketchsize64 == 0,
        "Sketch word count is not divisible by sketchsize64"
    );

    let bbits = num_bins / sketchsize64;

    // Same stride calculation as CUDA's flatten_by_samples().
    let sample_stride = 1usize;
    let bin_stride = n_samples.div_ceil(32) * 32;
    let kmer_stride = bin_stride
        .checked_mul(num_bins)
        .context("k-mer stride overflow")?;

    let total_words = kmer_stride
        .checked_mul(n_kmers)
        .context("Flattened sketch size overflow")?;

    // Zero-initialised, including the padded sample positions.
    let mut flattened = vec![0u64; total_words];

    // Prepare one slice per sample for the current k-mer length.
    // This avoids repeatedly calling get_sketch_slice() in the
    // innermost loop.
    for kmer_idx in 0..n_kmers {
        let slices: Vec<&[u64]> = (0..n_samples)
            .map(|sample_idx| {
                sketches.get_sketch_slice(sample_idx, kmer_idx)
            })
            .collect();

        for slice in &slices {
            ensure!(
                slice.len() == num_bins,
                "Inconsistent sketch word count"
            );
        }

        for bin_idx in 0..num_bins {
            let dst_base =
                kmer_idx * kmer_stride + bin_idx * bin_stride;

            for (sample_idx, slice) in slices.iter().enumerate() {
                flattened[dst_base + sample_idx] = slice[bin_idx];
            }
        }
    }

    // Convert to the u32 fields used by the WGSL uniform.
    let strides = SketchStrides {
        bin_stride: u32::try_from(bin_stride)?,
        kmer_stride: u32::try_from(kmer_stride)?,
        sample_stride: u32::try_from(sample_stride)?,
        sketchsize64: u32::try_from(sketchsize64)?,
        bbits: u32::try_from(bbits)?,
        _padding: [0; 3],
    };

    Ok((flattened, strides))
}

