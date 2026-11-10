// CUB host dispatch through a C ABI. All work uses caller-owned storage/stream.
#include <cuda_runtime.h>
#include <math_constants.h>
#include <cub/version.cuh>
#include <cub/device/device_batched_topk.cuh>
#include <cub/device/device_topk.cuh>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <cuda/argument>
#include <cuda/std/execution>
#include <cuda/stream_ref>
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>

static_assert(CUDA_VERSION >= 13050 && CUB_VERSION >= 300500,
              "GPU kNN requires CUDA Toolkit 13.5+ and CCCL/CUB 3.5+");

namespace {
constexpr size_t chunk_size = 4096;
constexpr size_t batched_k = 2048;
using Count = thrust::counting_iterator<size_t>;

size_t mul(size_t a, size_t b) {
    if (b && a > std::numeric_limits<size_t>::max() / b) throw std::overflow_error("size");
    return a * b;
}
size_t add(size_t a, size_t b) {
    if (a > std::numeric_limits<size_t>::max() - b) throw std::overflow_error("size");
    return a + b;
}
size_t aligned(size_t n) { return add(n, 255) & ~size_t(255); }
size_t reduced(size_t n, size_t k) {
    return add(mul(n / chunk_size, k), std::min(n % chunk_size, k));
}

// Modes: reference rows, query columns, condensed diagonal rows, intermediate rows.
struct View {
    const float* keys;
    const uint32_t* ids;
    size_t width, columns, samples;
    int mode;
    __host__ __device__ size_t neighbour(size_t row, size_t q) const {
        return mode == 2 && q >= row ? q + 1 : q;
    }
    __host__ __device__ size_t offset(size_t row, size_t q) const {
        if (mode == 0) return row * columns + q;
        if (mode == 1) return q * columns + row;
        if (mode == 3) return row * width + q;
        size_t j = neighbour(row, q);
        size_t a = row < j ? row : j, b = row < j ? j : row;
        return a * (2 * samples - a - 1) / 2 + b - a - 1;
    }
};
struct KeyAt {
    View view; size_t row, begin;
    __host__ __device__ float operator()(size_t q) const { return view.keys[view.offset(row, begin + q)]; }
};
struct IdAt {
    View view; size_t row, begin;
    __host__ __device__ uint32_t operator()(size_t q) const {
        return view.mode == 3 ? view.ids[view.offset(row, begin + q)]
                              : uint32_t(view.neighbour(row, begin + q));
    }
};
using Keys = thrust::transform_iterator<KeyAt, Count>;
using Ids = thrust::transform_iterator<IdAt, Count>;
struct KeyRows {
    View view; size_t chunks;
    __host__ __device__ Keys operator()(size_t segment) const {
        return Keys(Count(0), KeyAt{view, segment / chunks, (segment % chunks) * chunk_size});
    }
};
struct IdRows {
    View view; size_t chunks;
    __host__ __device__ Ids operator()(size_t segment) const {
        return Ids(Count(0), IdAt{view, segment / chunks, (segment % chunks) * chunk_size});
    }
};
template<class T> struct Outputs {
    T* ptr; size_t width, chunks, k;
    __host__ __device__ T* operator()(size_t segment) const {
        // Null pointers are used only by host-side workspace queries.
        return ptr ? ptr + (segment / chunks) * width + (segment % chunks) * k : nullptr;
    }
};
struct Lengths {
    size_t width, chunks;
    __host__ __device__ int operator()(size_t segment) const {
        size_t remaining = width - (segment % chunks) * chunk_size;
        return int(remaining < chunk_size ? remaining : chunk_size);
    }
};
auto environment(cudaStream_t stream) {
    return cuda::std::execution::env{cuda::stream_ref{stream}, cuda::execution::require(
        cuda::execution::determinism::not_guaranteed,
        cuda::execution::tie_break::unspecified,
        cuda::execution::output_ordering::unsorted)};
}
cudaError_t batch(void* temp, size_t& bytes, View view, size_t rows, size_t k,
                  float* out_keys, uint32_t* out_ids, cudaStream_t stream) {
    size_t chunks = (view.width - 1) / chunk_size + 1;
    size_t count = mul(rows, chunks);
    if (count > size_t(std::numeric_limits<int>::max())) return cudaErrorInvalidValue;
    size_t width = reduced(view.width, k);
    return cub::DeviceBatchedTopK::MinPairs(temp, bytes,
        thrust::make_transform_iterator(Count(0), KeyRows{view, chunks}),
        thrust::make_transform_iterator(Count(0), Outputs<float>{out_keys, width, chunks, k}),
        thrust::make_transform_iterator(Count(0), IdRows{view, chunks}),
        thrust::make_transform_iterator(Count(0), Outputs<uint32_t>{out_ids, width, chunks, k}),
        cuda::args::deferred_sequence{
            thrust::make_transform_iterator(Count(0), Lengths{view.width, chunks}),
            cuda::args::bounds<1, 4096>()},
        cuda::args::immediate{int(k), cuda::args::bounds<1, 2048>()},
        int(count), environment(stream));
}

cudaError_t direction(void* workspace, size_t& required, View initial, size_t rows,
                      size_t requested_k, float* out_keys, uint32_t* out_ids, cudaStream_t stream) {
    const size_t k = std::min(requested_k, initial.width);
    if (!rows || !k) { required = 0; return cudaSuccess; }
    size_t temp_bytes = 0;
    if (k > batched_k) {
        auto status = cub::DeviceTopK::MinPairs(nullptr, temp_bytes,
            Keys(Count(0), KeyAt{initial, 0, 0}), out_keys,
            Ids(Count(0), IdAt{initial, 0, 0}), out_ids,
            initial.width, k, environment(stream));
        if (status != cudaSuccess) return status;
        required = std::max(size_t(1), temp_bytes);
        if (!workspace) return cudaSuccess;
        for (size_t row = 0; row < rows; ++row) {
            size_t available = temp_bytes;
            status = cub::DeviceTopK::MinPairs(workspace, available,
                Keys(Count(0), KeyAt{initial, row, 0}), out_keys + row * k,
                Ids(Count(0), IdAt{initial, row, 0}), out_ids + row * k,
                initial.width, k, environment(stream));
            if (status != cudaSuccess) return status;
        }
        return cudaSuccess;
    }
    // Two ping-pong key/ID arrays accommodate every intermediate reduction.
    const size_t buffer_bytes = initial.width > chunk_size
        ? aligned(mul(mul(rows, reduced(initial.width, k)), sizeof(float))) : 0;
    const size_t prefix = mul(4, buffer_bytes);
    View query = initial;
    for (;;) {
        size_t bytes = 0;
        auto status = batch(nullptr, bytes, query, rows, k, nullptr, nullptr, stream);
        if (status != cudaSuccess) return status;
        temp_bytes = std::max(temp_bytes, bytes);
        if (query.width <= chunk_size) break;
        query.width = reduced(query.width, k);
        query.mode = 3;
    }
    required = add(prefix, std::max(size_t(1), temp_bytes));
    if (!workspace) return cudaSuccess;
    auto* base = static_cast<unsigned char*>(workspace);
    View view = initial;
    size_t iteration = 0;
    for (;;) {
        bool final = view.width <= chunk_size;
        size_t slot = iteration++ % 2;
        float* keys = final ? out_keys : reinterpret_cast<float*>(base + slot * 2 * buffer_bytes);
        uint32_t* ids = final ? out_ids : reinterpret_cast<uint32_t*>(base + (slot * 2 + 1) * buffer_bytes);
        size_t available = temp_bytes;
        auto status = batch(base + prefix, available, view, rows, k, keys, ids, stream);
        if (status != cudaSuccess) return status;
        if (final) return cudaSuccess;
        view = View{keys, ids, reduced(view.width, k), 0, 0, 3};
    }
}
} // namespace

extern "C" int sketchlib_tile_topk(void* workspace, size_t* bytes, const float* distances,
    float* out_keys, uint32_t* out_ids, size_t nr, size_t nq, size_t k, int diagonal, void* raw_stream) noexcept {
    try {
        if (!bytes || !nr || !nq || !k || nr > UINT32_MAX || nq > UINT32_MAX
            || (diagonal && nr != nq)) return cudaErrorInvalidValue;
        auto stream = static_cast<cudaStream_t>(raw_stream);
        View refs{distances, nullptr, diagonal ? nr - 1 : nq, nq, nr, diagonal ? 2 : 0};
        View queries{distances, nullptr, nr, nq, nr, 1};
        size_t ref_bytes = 0, query_bytes = 0;
        auto status = direction(nullptr, ref_bytes, refs, nr, k, nullptr, nullptr, stream);
        if (status != cudaSuccess) return status;
        if (!diagonal) {
            status = direction(nullptr, query_bytes, queries, nq, k, nullptr, nullptr, stream);
            if (status != cudaSuccess) return status;
        }
        size_t needed = std::max(size_t(1), std::max(ref_bytes, query_bytes));
        if (!workspace) { *bytes = needed; return cudaSuccess; }
        if (*bytes < needed || !distances || !out_keys || !out_ids) return cudaErrorInvalidValue;
        status = direction(workspace, ref_bytes, refs, nr, k, out_keys, out_ids, stream);
        if (status != cudaSuccess || diagonal) return status;
        size_t offset = mul(nr, std::min(k, nq));
        return direction(workspace, query_bytes, queries, nq, k, out_keys + offset, out_ids + offset, stream);
    } catch (const std::overflow_error&) { return cudaErrorInvalidValue; }
      catch (...) { return cudaErrorUnknown; }
}

extern "C" const char* sketchlib_topk_error_string(int status) noexcept {
    return cudaGetErrorString(static_cast<cudaError_t>(status));
}

namespace {
// A row reserves min(k, block_size) entries for each neighbour block. Its own
// block contributes one fewer entry when block_size <= k; only that slot needs
// initialisation. All other entries are written exactly once by tile scatters.
__device__ size_t device_min(size_t a, size_t b) { return a < b ? a : b; }
__global__ void pad_candidates(float* keys, uint32_t* ids, size_t n, size_t block,
                               size_t k, size_t width) {
    size_t row = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    if (row >= n) return;
    size_t owner = row / block;
    size_t samples = device_min(block, n - owner * block);
    if (samples <= k) {
        size_t offset = row * width + owner * device_min(k, block) + samples - 1;
        keys[offset] = CUDART_INF_F;
        ids[offset] = UINT32_MAX;
    }
}

__global__ void scatter_candidates(const float* keys, const uint32_t* ids,
    float* all_keys, uint32_t* all_ids, size_t nr, size_t ref_k, size_t query_k,
    size_t ref_start, size_t query_start, size_t block, size_t k, size_t width,
    size_t count) {
    size_t i = size_t(blockIdx.x) * blockDim.x + threadIdx.x;
    // Grid-stride iteration also handles candidate counts larger than a grid.
    for (; i < count; i += size_t(blockDim.x) * gridDim.x) {
        size_t refs = nr * ref_k;
        bool reference = i < refs;
        size_t local = reference ? i : i - refs;
        size_t row_k = reference ? ref_k : query_k;
        size_t start = reference ? ref_start : query_start;
        size_t neighbour = reference ? query_start : ref_start;
        size_t offset = (start + local / row_k) * width
            + (neighbour / block) * device_min(k, block) + local % row_k;
        all_keys[offset] = keys[i];
        all_ids[offset] = uint32_t(neighbour + ids[i]);
    }
}
} // namespace

extern "C" int sketchlib_init_candidates(float* keys, uint32_t* ids, size_t n,
    size_t block, size_t k, size_t width, void* raw_stream) noexcept {
    if (!keys || !ids || !n || !block || !k || !width) return cudaErrorInvalidValue;
    pad_candidates<<<(n + 255) / 256, 256, 0, static_cast<cudaStream_t>(raw_stream)>>>(
        keys, ids, n, block, k, width);
    return cudaGetLastError();
}

extern "C" int sketchlib_scatter_candidates(const float* keys, const uint32_t* ids,
    float* all_keys, uint32_t* all_ids, size_t nr, size_t nq, size_t k,
    size_t ref_start, size_t query_start, size_t block, size_t width,
    void* raw_stream) noexcept {
    try {
        if (!keys || !ids || !all_keys || !all_ids || !nr || !nq || !k || !block || !width)
            return cudaErrorInvalidValue;
        bool diagonal = ref_start == query_start;
        size_t ref_k = std::min(k, diagonal ? nr - 1 : nq);
        size_t query_k = diagonal ? 0 : std::min(k, nr);
        size_t count = add(mul(nr, ref_k), mul(nq, query_k));
        if (!count) return cudaSuccess;
        unsigned grid = unsigned(std::min(size_t(65535), (count - 1) / 256 + 1));
        scatter_candidates<<<grid, 256, 0, static_cast<cudaStream_t>(raw_stream)>>>(
            keys, ids, all_keys, all_ids, nr, ref_k, query_k, ref_start, query_start,
            block, k, width, count);
        return cudaGetLastError();
    } catch (...) { return cudaErrorInvalidValue; }
}

extern "C" int sketchlib_global_topk(void* workspace, size_t* bytes,
    const float* keys, const uint32_t* ids, float* out_keys, uint32_t* out_ids,
    size_t rows, size_t width, size_t k, void* raw_stream) noexcept {
    try {
        if (!bytes || !rows || !k || k > width) return cudaErrorInvalidValue;
        mul(rows, width);
        mul(rows, k);
        View view{keys, ids, width, 0, 0, 3};
        size_t needed = 0;
        auto stream = static_cast<cudaStream_t>(raw_stream);
        auto status = direction(nullptr, needed, view, rows, k, nullptr, nullptr, stream);
        if (status != cudaSuccess) return status;
        if (!workspace) { *bytes = needed; return cudaSuccess; }
        if (*bytes < needed || !keys || !ids || !out_keys || !out_ids)
            return cudaErrorInvalidValue;
        return direction(workspace, needed, view, rows, k, out_keys, out_ids, stream);
    } catch (const std::overflow_error&) { return cudaErrorInvalidValue; }
      catch (...) { return cudaErrorUnknown; }
}
