// Host [sample][k-mer][bin] -> device [k-mer][bin][sample], for one tile range, as the
// CUDA code's interleave_sketches(). Samples are padded with zeros up to bin_stride (a
// multiple of 32, as pp-sketchlib's warpPad())
struct InterleaveParams {
    n_samples: u32,   // samples in this tile range
    n_kmers: u32,
    n_bins: u32,      // words per k-mer sketch, all bit planes included
    bin_stride: u32,  // the tile-range size, padded to a multiple of 32
};

@group(0) @binding(0) var<uniform> params: InterleaveParams;
@group(0) @binding(1) var<storage, read> sketch_slice: array<vec2<u32>>;
@group(0) @binding(2) var<storage, read_write> interleaved: array<vec2<u32>>;

override WORKGROUP_SIZE: u32 = 256;

@compute @workgroup_size(WORKGROUP_SIZE, 1, 1)
fn main(
    @builtin(global_invocation_id) global_invocation_id: vec3<u32>,
    @builtin(num_workgroups) num_workgroups: vec3<u32>,
) {
    let total = params.n_kmers * params.n_bins * params.bin_stride;
    let sample_stride_in = params.n_bins * params.n_kmers;
    // Grid-stride loop (as pp-sketchlib's dist.cu:141-142), so one dispatch covers any size
    for (var idx = global_invocation_id.x; idx < total; idx += num_workgroups.x * WORKGROUP_SIZE) {
        // In the output, k-mer outer index, then bin, then sample
        let sample_idx = idx % params.bin_stride;
        let bin_idx = (idx / params.bin_stride) % params.n_bins;
        let kmer_idx = idx / (params.bin_stride * params.n_bins);

        var word = vec2<u32>(0u);
        if (sample_idx < params.n_samples) {
            word = sketch_slice[sample_idx * sample_stride_in + kmer_idx * params.n_bins + bin_idx];
        }
        interleaved[idx] = word;
    }
}
