// CUB host dispatch through a C ABI. All work uses caller-owned storage/stream.
#include <cuda_runtime.h>
#include <cub/version.cuh>
#include <cub/device/device_batched_topk.cuh>
#include <cub/device/device_topk.cuh>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <cuda/argument>
#include <cuda/std/execution>
#include <cuda/stream_ref>
#include <cstdint>

static_assert(CUDA_VERSION >= 13050 && CUB_VERSION >= 300500,
              "GPU kNN requires CUDA Toolkit 13.5+ and CCCL/CUB 3.5+");

namespace {
// Wrapper policy, not universal CUB limits: CCCL 3.5's batched implementation
// requires each segment to fit one thread block's resources. These bounds and
// the annotations in the batched adapter select that implementation for our float/u32 pairs.
constexpr size_t chunk_size = 8192;
constexpr size_t batched_k = 2048;
enum class Layout : int {
    ReferenceRows = 0,
    QueryColumns = 1,
    DiagonalRows = 2,
};
using Count = thrust::counting_iterator<size_t>;

// Logical rows over existing storage; no transposed/full diagonal matrix is made.
// width excludes the self-pair for DiagonalRows;
// columns is the original rectangular stride, samples the diagonal sample count.
struct View {
    const float* keys;
    size_t width, columns, samples;
    Layout layout;
    __host__ __device__ size_t neighbour(size_t row, size_t q) const {
        // Expand the self-excluding logical position to its block-local ID.
        return layout == Layout::DiagonalRows && q >= row ? q + 1 : q;
    }
    __host__ __device__ size_t offset(size_t row, size_t q) const {
        if (layout == Layout::ReferenceRows) return row * columns + q;
        if (layout == Layout::QueryColumns) return q * columns + row;
        size_t j = neighbour(row, q);
        size_t a = row < j ? row : j;
        size_t b = row < j ? j : row;
        return a * (2 * samples - a - 1) / 2 + b - a - 1;
    }
};
struct KeyAt {
    View view; size_t row, begin;
    __host__ __device__ float operator()(size_t q) const {
        return view.keys[view.offset(row, begin + q)];
    }
};
struct IdAt {
    View view; size_t row, begin;
    __host__ __device__ uint32_t operator()(size_t q) const {
        return uint32_t(view.neighbour(row, begin + q));
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
    // Selection order and which equal-valued cutoff candidates survive are
    // unspecified. The host heaps perform the final global selection/order.
    return cuda::std::execution::env{cuda::stream_ref{stream}, cuda::execution::require(
        cuda::execution::determinism::not_guaranteed,
        cuda::execution::tie_break::unspecified,
        cuda::execution::output_ordering::unsorted)};
}
} // namespace

// Trusted Rust callers validate dimensions and bind the owning device. Null
// workspace queries CUB scratch; execution uses the supplied capacity directly.
// All pointers/stream remain caller-owned, with no allocation or synchronization.
extern "C" void sketchlib_topk_limits(size_t* chunk, size_t* k) noexcept {
    *chunk = chunk_size;
    *k = batched_k;
}

extern "C" int sketchlib_batched_topk(void* temp, size_t* bytes, const View* input,
    size_t chunks, size_t output_width, int k, int segments,
    float* out_keys, uint32_t* out_ids, void* raw_stream) noexcept {
    try {
        View view = *input;
        return cub::DeviceBatchedTopK::MinPairs(temp, *bytes,
            thrust::make_transform_iterator(Count(0), KeyRows{view, chunks}),
            thrust::make_transform_iterator(Count(0), Outputs<float>{out_keys, output_width, chunks, size_t(k)}),
            thrust::make_transform_iterator(Count(0), IdRows{view, chunks}),
            thrust::make_transform_iterator(Count(0), Outputs<uint32_t>{out_ids, output_width, chunks, size_t(k)}),
            cuda::args::deferred_sequence{
                thrust::make_transform_iterator(Count(0), Lengths{view.width, chunks}),
                cuda::args::bounds<1, int(chunk_size)>()},
            cuda::args::immediate{k, cuda::args::bounds<1, int(batched_k)>()},
            segments, environment(static_cast<cudaStream_t>(raw_stream)));
    } catch (...) {
        return cudaErrorUnknown;
    }
}

extern "C" int sketchlib_rows_topk(void* temp, size_t* bytes, const View* input,
    size_t rows, size_t k, float* out_keys, uint32_t* out_ids, void* raw_stream) noexcept {
    try {
        View view = *input;
        auto stream = static_cast<cudaStream_t>(raw_stream);
        if (!temp) {
            return cub::DeviceTopK::MinPairs(nullptr, *bytes,
                Keys(Count(0), KeyAt{view, 0, 0}), out_keys,
                Ids(Count(0), IdAt{view, 0, 0}), out_ids,
                view.width, k, environment(stream));
        }
        // Keep this short loop beside CUB to avoid one FFI call per row. Each row
        // reuses caller-provided scratch, ordered by the same stream.
        for (size_t row = 0; row < rows; ++row) {
            size_t available = *bytes;
            auto status = cub::DeviceTopK::MinPairs(temp, available,
                Keys(Count(0), KeyAt{view, row, 0}), out_keys + row * k,
                Ids(Count(0), IdAt{view, row, 0}), out_ids + row * k,
                view.width, k, environment(stream));
            if (status != cudaSuccess) return status;
        }
        return cudaSuccess;
    } catch (...) {
        return cudaErrorUnknown;
    }
}

extern "C" const char* sketchlib_topk_error_string(int status) noexcept {
    return cudaGetErrorString(static_cast<cudaError_t>(status));
}
