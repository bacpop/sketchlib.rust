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
#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>

static_assert(CUDA_VERSION >= 13050 && CUB_VERSION >= 300500,
              "GPU kNN requires CUDA Toolkit 13.5+ and CCCL/CUB 3.5+");

namespace {
// Wrapper policy, not universal CUB limits: CCCL 3.5's batched implementation
// requires each segment to fit one thread block's resources. These bounds and
// the annotations in batch() select that implementation for our float/u32 pairs.
constexpr size_t chunk_size = 4096;
constexpr size_t batched_k = 2048;
constexpr size_t workspace_alignment = 256;
static_assert(sizeof(float) == sizeof(uint32_t),
              "Intermediate key and ID buffers must have equal element sizes");

enum class Layout : int {
    ReferenceRows = 0,
    QueryColumns = 1,
    DiagonalRows = 2,
    IntermediateRows = 3,
};
using Count = thrust::counting_iterator<size_t>;

size_t mul(size_t a, size_t b) {
    if (b && a > std::numeric_limits<size_t>::max() / b) {
        throw std::overflow_error("size");
    }
    return a * b;
}
size_t add(size_t a, size_t b) {
    if (a > std::numeric_limits<size_t>::max() - b) {
        throw std::overflow_error("size");
    }
    return a + b;
}
size_t aligned(size_t n) {
    return add(n, workspace_alignment - 1) & ~(workspace_alignment - 1);
}
size_t reduced(size_t n, size_t k) {
    return add(mul(n / chunk_size, k), std::min(n % chunk_size, k));
}

// Logical rows over existing storage; no transposed/full diagonal matrix is made.
// width excludes the self-pair for DiagonalRows;
// columns is the original rectangular stride, samples the diagonal sample count.
struct View {
    const float* keys;
    const uint32_t* ids;
    size_t width, columns, samples;
    Layout layout;
    __host__ __device__ size_t neighbour(size_t row, size_t q) const {
        // Expand the self-excluding logical position to its block-local ID.
        return layout == Layout::DiagonalRows && q >= row ? q + 1 : q;
    }
    __host__ __device__ size_t offset(size_t row, size_t q) const {
        if (layout == Layout::ReferenceRows) return row * columns + q;
        if (layout == Layout::QueryColumns) return q * columns + row;
        if (layout == Layout::IntermediateRows) return row * width + q;
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
        // Reductions carry original neighbour IDs, not intermediate positions.
        return view.layout == Layout::IntermediateRows
            ? view.ids[view.offset(row, begin + q)]
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
    // Selection order and which equal-valued cutoff candidates survive are
    // unspecified. The host heaps perform the final global selection/order.
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
            cuda::args::bounds<1, int(chunk_size)>()},
        cuda::args::immediate{int(k), cuda::args::bounds<1, int(batched_k)>()},
        int(count), environment(stream));
}

cudaError_t direction(void* workspace, size_t& required, View initial, size_t rows,
                      size_t requested_k, float* out_keys, uint32_t* out_ids, cudaStream_t stream) {
    const size_t k = std::min(requested_k, initial.width);
    if (!rows || !k) {
        required = 0;
        return cudaSuccess;
    }
    size_t temp_bytes = 0;
    // Larger effective k bypasses our batched bound. Rows reuse one DeviceTopK
    // scratch allocation sequentially on the caller's stream.
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
    // Keeping each chunk's smallest k preserves a valid smallest k for the
    // whole row, including arbitrary choices among equal-valued cutoff ties.
    // Repeat until one batched segment can select the row's final candidates.
    //
    // Workspace prefix: keys0, ids0, keys1, ids1, each buffer_bytes bytes and
    // 256-byte aligned (float and uint32_t are both four bytes). The first
    // reduction is the largest; later passes alternate these two buffer pairs.
    // CUB scratch follows the prefix and is sized for the largest pass request.
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
        query.layout = Layout::IntermediateRows;
    }
    required = add(prefix, std::max(size_t(1), temp_bytes));
    if (!workspace) return cudaSuccess;
    auto* base = static_cast<unsigned char*>(workspace);
    View view = initial;
    size_t iteration = 0;
    for (;;) {
        bool final = view.width <= chunk_size;
        size_t slot = iteration++ % 2;
        // Final candidates go directly to the caller's output, not the prefix.
        float* keys = final
            ? out_keys
            : reinterpret_cast<float*>(base + slot * 2 * buffer_bytes);
        uint32_t* ids = final
            ? out_ids
            : reinterpret_cast<uint32_t*>(base + (slot * 2 + 1) * buffer_bytes);
        size_t available = temp_bytes;
        auto status = batch(base + prefix, available, view, rows, k, keys, ids, stream);
        if (status != cudaSuccess) return status;
        if (final) return cudaSuccess;
        view = View{keys, ids, reduced(view.width, k), 0, 0, Layout::IntermediateRows};
    }
}
} // namespace

// The caller must bind the owning device before either mode. A null workspace
// queries required bytes without dereferencing device pointers or enqueuing work.
// Otherwise *bytes is the available workspace capacity; all storage belongs to
// the caller and must remain valid until the supplied stream completes. This
// wrapper neither allocates device memory nor synchronizes the stream.
extern "C" int sketchlib_tile_topk(void* workspace, size_t* bytes, const float* distances,
    float* out_keys, uint32_t* out_ids, size_t nr, size_t nq, size_t k, int diagonal, void* raw_stream) noexcept {
    try {
        if (!bytes || !nr || !nq || !k || nr > UINT32_MAX || nq > UINT32_MAX
            || (diagonal && nr != nq)) {
            return cudaErrorInvalidValue;
        }
        auto stream = static_cast<cudaStream_t>(raw_stream);
        View refs{distances, nullptr, diagonal ? nr - 1 : nq, nq, nr,
                  diagonal ? Layout::DiagonalRows : Layout::ReferenceRows};
        View queries{distances, nullptr, nr, nq, nr, Layout::QueryColumns};
        size_t ref_bytes = 0, query_bytes = 0;
        auto status = direction(nullptr, ref_bytes, refs, nr, k, nullptr, nullptr, stream);
        if (status != cudaSuccess) return status;
        if (!diagonal) {
            status = direction(nullptr, query_bytes, queries, nq, k, nullptr, nullptr, stream);
            if (status != cudaSuccess) return status;
        }
        size_t needed = std::max(size_t(1), std::max(ref_bytes, query_bytes));
        if (!workspace) {
            *bytes = needed;
            return cudaSuccess;
        }
        if (*bytes < needed || !distances || !out_keys || !out_ids) {
            return cudaErrorInvalidValue;
        }
        // Both directions run on the same stream, so workspace is their maximum,
        // not their sum. Outputs are unsorted rows of min(k, available neighbours):
        // nr reference rows first, then nq query rows; diagonals need only refs.
        status = direction(workspace, ref_bytes, refs, nr, k, out_keys, out_ids, stream);
        if (status != cudaSuccess || diagonal) return status;
        size_t offset = mul(nr, std::min(k, nq));
        return direction(workspace, query_bytes, queries, nq, k, out_keys + offset, out_ids + offset, stream);
    } catch (const std::overflow_error&) {
        return cudaErrorInvalidValue;
    } catch (...) {
        return cudaErrorUnknown;
    }
}

extern "C" const char* sketchlib_topk_error_string(int status) noexcept {
    return cudaGetErrorString(static_cast<cudaError_t>(status));
}
