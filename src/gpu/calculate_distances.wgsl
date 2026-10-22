// Initial struct declarations
struct SketchStrides {
    bin_stride: u32,     // 4 bytes, offset 0
    kmer_stride: u32,    // 4 bytes, offset 4
    sample_stride: u32,  // 4 bytes, offset 8
    sketchsize64: u32,   // 4 bytes, offset 12
    bbits: u32,          // 4 bytes, offset 16
    _pad0: u32,          // 4 bytes, offset 20
    _pad1: u32,          // 4 bytes, offset 24
    _pad2: u32,          // 4 bytes, offset 28
};

struct DistanceParams {
    ref_strides: SketchStrides,   // 32 bytes, offset 0
    query_strides: SketchStrides, // 32 bytes, offset 32

    ref_n: u32,           // 4 bytes, offset 64
    query_n: u32,         // 4 bytes, offset 68
    kmer_n: u32,          // 4 bytes, offset 72
    dist_n: u32,          // 4 bytes, offset 76
    self_comparison: u32, // 4 bytes, offset 80
    use_shared: u32,      // 4 bytes, offset 84

    _pad1: u32,          // 4 bytes, offset 88
    _pad2: u32,          // 4 bytes, offset 92
};

// Get your bindings
@group(0) @binding(0) var<uniform> params: DistanceParams;
@group(0) @binding(1) var<storage, read> input_sketches: array<vec2<u32>>;
// This storage buffer can be read from and written to
@group(0) @binding(2) var<storage, read_write> output_dists: array<f32>;
// k-mer lengths, CUDA's `const int *kmers`
@group(0) @binding(3) var<storage, read> kmers: array<u32>;
// Completeness of each sample (only read when USE_COMPLETENESS)
@group(0) @binding(4) var<storage, read> completeness: array<f32>;

// Tells wgpu that this function is a valid compute pipeline entry_point
override WORKGROUP_SIZE: u32 = 256;
// Words of the query sketch held in workgroup memory (1 when it doesn't fit): CUDA's
// shared memory size at launch
override SHARED_WORDS: u32 = 1u;
// Invocations that copy the query into workgroup memory: CUDA's warp_size
override QUERY_LOADERS: u32 = 32u;
// Not in CUDA: CPU-only features
override ANI: bool = false;
override USE_COMPLETENESS: bool = false;
override COMPLETENESS_CUTOFF: f32 = 0.64;

// CUDA's extern __shared__ uint64_t query_shared[]
var<workgroup> query_shared: array<vec2<u32>, SHARED_WORDS>;


// Specifies the "dimension" of this work group
// It cannot be accessed later! Thus, I define a constant that I can use in the different functions
// The compute attribute defines that this is a compute shader. IF YOU PUT THINGS BETWEEN IT AND @WORKGROUP_SIZE, BAD THINGS HAPPEN!!!
@compute @workgroup_size(WORKGROUP_SIZE, 1, 1)

// NOTE: because I keep forgetting
// CUDA ================ WGSL
// thread -------------- (local) invocation
// warp   -------------- subgroup
// block  -------------- workgroup
// grid   -------------- workgroups
// global thread index - global invocation index


fn main(
    // global_invocation_id specifies our position in the invocation grid
    // builtin is the way of importting special values from the execution environment, normally these workgroup/invocation IDs
    @builtin(global_invocation_id) global_invocation_id: vec3<u32>,
    @builtin(workgroup_id) workgroup_id: vec3<u32>,
    @builtin(local_invocation_id) local_invocation_id: vec3<u32>,
) {
    // This follows the CUDA kernel indexing FOR QUERY VS REF
//     let query_idx   = workgroup_id.x / BLOCKS_PER_QUERY;
//     let ref_idx     = (workgroup_id.x % BLOCKS_PER_QUERY) * WORKGROUP_SIZE + local_invocation_id.x;
//     let dist_idx    = query_idx * REF_N + ref_idx;

    // This is for ALL-vs-ALL. This gets "query_idx", which allows you to know which query genome we are processing
    var workgroups_done : u32 = 0;
    var query_idx: u32 = 0;
    loop {
        if (query_idx >= params.ref_n - 1u) { // Suggestion from CHatty to prevent underflows in veeeeeeeery large data sets (should not happen ever... I think?)
            break;
        }

        workgroups_done +=
            (params.ref_n + WORKGROUP_SIZE - 2 - query_idx) / WORKGROUP_SIZE;

        if (workgroups_done > workgroup_id.x) {
            break;
        }
        query_idx++;
    }

    // How many workgroups/blocks are there per query (i.e. share the same query genome in their shared memory)
    let workgroups_per_query : u32 = (params.ref_n + WORKGROUP_SIZE - 2 - query_idx) / WORKGROUP_SIZE;
    // Local ID in the query, starting from zero
    let local_workgroup_idx : u32 = workgroup_id.x - (workgroups_done - workgroups_per_query);
    // What reference genome we are going to calculate the distance against
    let ref_idx : u32 = query_idx + 1u + local_workgroup_idx * WORKGROUP_SIZE + local_invocation_id.x;

    // SO that we don't access references that are not there! Invocations past the last
    // reference must still reach every workgroupBarrier(), so no early return
    let valid = ref_idx < params.ref_n;
    let dist_idx =
        query_idx * (2u * params.ref_n - query_idx - 1u) / 2u
        + (ref_idx - query_idx - 1u);

    // Entries of the flatten vector where we have the data
    // (the binding already starts at the selected k-mer, as CUDA's `ref + k * kmer_stride`)
    var ref_start = ref_idx * params.ref_strides.sample_stride;
    var query_start = query_idx * params.query_strides.sample_stride;
    let query_words = params.query_strides.sketchsize64 * params.query_strides.bbits;

    // Loop over k-mer lengths, as CUDA's calculate_dists. Jaccard runs a single k-mer
    var dist: f32 = 0.0;
    for (var kmer_idx = 0u; kmer_idx < params.kmer_n; kmer_idx++) {
        // Nobody may still be reading the previous k-mer's query (CUDA's block.sync())
        workgroupBarrier();
        // Copy the query sketch into workgroup memory with QUERY_LOADERS invocations
        if (params.use_shared != 0u && local_invocation_id.x < QUERY_LOADERS) {
            for (var w = local_invocation_id.x; w < query_words; w += QUERY_LOADERS) {
                query_shared[w] =
                    input_sketches[query_start + w * params.query_strides.bin_stride];
            }
        }
        // CUDA's barrier.arrive_and_wait()
        workgroupBarrier();

        if (valid) {
            let samebits = count_samebits(ref_start, query_start);
            dist = jaccard_output(samebits, query_idx, ref_idx, kmers[kmer_idx]);
        }

        // Move to next k-mer length
        ref_start += params.ref_strides.kmer_stride;
        query_start += params.query_strides.kmer_stride;
    }

    if (valid) {
        output_dists[dist_idx] = dist;
    }
}

// A query word, from workgroup memory or from the storage buffer. CUDA switches
// query_ptr between the two; WGSL can't point at either, so branch on the flag
fn query_word(query_start: u32, bin_index: u32) -> vec2<u32> {
    if (params.use_shared != 0u) {
        return query_shared[bin_index];
    }
    return input_sketches[query_start + bin_index * params.query_strides.bin_stride];
}

// Number of matching bins (CUDA's jaccard_dist() up to the popcount)
fn count_samebits(ref_start: u32, query_start: u32) -> u32 {
    var samebits = 0u;

    for (var i = 0u; i < params.ref_strides.sketchsize64; i++) {
        var bin_index = i * params.ref_strides.bbits;
        var bits = vec2<u32>(0xffffffffu); // all 1s

        for (var j = 0u; j < params.ref_strides.bbits; j++) {
            let ref_word =
                input_sketches[ref_start + bin_index * params.ref_strides.bin_stride];

            bits &= ~(ref_word ^ query_word(query_start, bin_index));
            bin_index++;
        }

        samebits +=
            countOneBits(bits[0]) + countOneBits(bits[1]);
    }
    return samebits;
}

// jaccard_index_generic() plus the distance/ANI step of the CPU code
fn jaccard_output(samebits: u32, query_idx: u32, ref_idx: u32, kmer: u32) -> f32 {
    // Correction for randomness
    let total_bins = f32(params.ref_strides.sketchsize64 * 64u);
    let expected_random = total_bins / f32(1u << params.ref_strides.bbits);
    let denominator = total_bins - expected_random; // exact in f32
    var jaccard = clamp((f32(samebits) - expected_random) / denominator, 0.0, 1.0);

    var corrected = false;
    if (USE_COMPLETENESS) {
        let c1 = completeness[query_idx];
        let c2 = completeness[ref_idx];
        if (c1 * c2 >= COMPLETENESS_CUTOFF) {
            jaccard = min(jaccard / (c1 * c2 / (c1 + c2 - c1 * c2)), 1.0);
            corrected = true;
        }
    }

    if (ANI) {
        // log(0) is undefined in WGSL; on the CPU ln(0) = -inf ends as max(0, -inf) = 0
        if (jaccard == 0.0) {
            return 0.0;
        }
        return max(0.0, 1.0 + 1.0 / f32(kmer) * log(2.0 * jaccard / (1.0 + jaccard)));
    }
    if (corrected) {
        return 1.0 - jaccard;
    }
    // 1 - J without cancellation: total_bins - samebits is exact
    return clamp((total_bins - f32(samebits)) / denominator, 0.0, 1.0);
}
