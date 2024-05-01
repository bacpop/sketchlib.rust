/*
 *
 * dist.cu
 * Sketchlib dists using CUDA
 *
 */

// std
#include <algorithm>
#include <assert.h>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <tuple>
#include <unistd.h>
#include <vector>

// memcpy_async
#include <cub/cub.cuh>
#include <cuda/barrier>
#include <cooperative_groups.h>
#pragma nv_diag_suppress static_var_with_dynamic_init

// internal headers
#include "cuda.cuh"
#include "containers.cuh"
#include "matrix_idx.hpp"
#include "gpu.hpp"

#define BITATPOS(x, pos) ((x & (1ULL << pos)) >> pos)
#define NBITS(x) (8 * sizeof(x))
#define ROUNDDIV(a, b) (((a) + (b) / 2) / (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/******************
 *			          *
 *	Device code   *
 *			          *
 *******************/

// Ternary used in observed_excess
template <class T> __device__ T non_neg_minus(T a, T b) {
  return a > b ? (a - b) : 0;
}

// Calculates excess observations above a random value
template <class T> __device__ T observed_excess(T obs, T exp, T max) {
  T diff = non_neg_minus(obs, exp);
  return (diff * max / (max - exp));
}

// CUDA version of bindash dist function (see dist.cpp)
__device__ float jaccard_dist(const uint64_t *sketch1, const uint64_t *sketch2,
                              const size_t sketchsize64, const size_t bbits,
                              const size_t s1_stride, const size_t s2_stride) {
  size_t samebits = 0;
  for (int i = 0; i < sketchsize64; i++) {
    int bin_index = i * bbits;
    uint64_t bits = ~((uint64_t)0ULL);
    for (int j = 0; j < bbits; j++) {
      // Almost all kernel time is spent on this line
      // (bbits * sketchsize64 * N^2 * 2 8-byte memory loads)
      bits &=
          ~(sketch1[bin_index * s1_stride] ^ sketch2[bin_index * s2_stride]);
      bin_index++;
    }
    samebits += __popcll(bits); // CUDA 64-bit popcnt
  }

  const size_t maxnbits = sketchsize64 * NBITS(uint64_t);
  const size_t expected_samebits = (maxnbits >> bbits);
  size_t intersize = samebits;
  if (!expected_samebits) {
    size_t ret = observed_excess(samebits, expected_samebits, maxnbits);
  }

  size_t unionsize = NBITS(uint64_t) * sketchsize64;
  float jaccard = __fdiv_ru(intersize, unionsize);
  return (jaccard);
}

// Simple linear regression, exact solution
__device__ void simple_linear_regression(float dists[],
                                         const float xsum, const float ysum,
                                         const float xysum,
                                         const float xsquaresum,
                                         const float ysquaresum, const int n,
                                         const int dist_col) {
  if (n < 2) {
    dists[0] = 0.0f;
    dists[1] = 0.0f;
  } else {
    // CUDA fast-math intrinsics on floats, which give comparable accuracy
    // Speed gain is fairly minimal, as most time spent on Jaccard distance
    // __fmul_ru(x, y) = x * y and rounds up.
    // __fpow(x, a) = x^a give 0 for x<0, so not using here (and it is slow)
    float xbar = xsum / n;
    float ybar = ysum / n;
    float x_diff = xsquaresum - __fmul_ru(xsum, xsum) / n;
    float y_diff = ysquaresum - __fmul_ru(ysum, ysum) / n;
    float xstddev = __fsqrt_ru((xsquaresum - __fmul_ru(xsum, xsum) / n) / n);
    float ystddev = __fsqrt_ru((ysquaresum - __fmul_ru(ysum, ysum) / n) / n);
    float r =
        __fdiv_ru(xysum - __fmul_ru(xsum, ysum) / n, __fsqrt_ru(x_diff * y_diff));
    float beta = __fmul_ru(r, __fdiv_ru(ystddev, xstddev));
    float alpha = __fmaf_ru(-beta, xbar, ybar); // maf: x * y + z

    // Store core/accessory in dists, truncating at zero
    // Memory should be initialised to zero so else block not strictly
    // necessary, but better safe than sorry!
    if (beta < 0.0f) {
      dists[0] = 1.0f - __expf(beta);
    } else {
      dists[0] = 0.0f;
    }

    if (alpha < 0.0f) {
      dists[1] = 1.0f - __expf(alpha);
    } else {
      dists[1] = 0.0f;
    }
  }
}

/******************
 *			      *
 *	Global code   *
 *			      *
 *******************/

// Main kernel functions run on the device,
// but callable from the host

__global__ void set_idx(long* idx, size_t row_samples, size_t col_samples, size_t col_offset) {
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < row_samples * col_samples;
    i += blockDim.x * gridDim.x) {
    idx[i] = col_offset + i % col_samples;
  }
}

__global__ void copy_top_k(float* sorted_dists, long* sorted_idx,
  float* all_sorted_dists, long* all_sorted_idx, int segment_size, int n_out,
    int kNN, int sketch_block_idx, int n_chunks, bool second_sort) {
  for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < n_out;
    i += blockDim.x * gridDim.x) {
    const int offset_in = (i / kNN) * segment_size + i % kNN;
    int offset_out = i;
    // If copying from the sorted kNN * n_chunk list into the final sparse matrix
    if (!second_sort) {
      // If copying from the sorted n_chunk dense list into the kNN * n_chunk list
      offset_out += ((i / kNN) * (kNN * (n_chunks - 1))) + (kNN * sketch_block_idx);
    }
    all_sorted_dists[offset_out] = sorted_dists[offset_in];
    all_sorted_idx[offset_out] = sorted_idx[offset_in];
  }
}

__global__ void calculate_dists(
    const bool self, const uint64_t *ref, const long ref_n,
    const uint64_t *query, const long query_n, const int *kmers,
    const int kmer_n, float *dists, const long long dist_n,
    const float *random_table, const uint16_t *ref_idx_lookup,
    const uint16_t *query_idx_lookup, const SketchStrides ref_strides,
    const SketchStrides query_strides, const RandomStrides random_strides,
    progress_ptrs progress, const bool use_shared,
    const int dist_col, const bool max_diagonal) {
  // Calculate indices for query, ref and results
  int ref_idx, query_idx, dist_idx;
  if (self) {
    // Blocks have the same i -- calculate blocks needed by each row up
    // to this point (blockIdx.x)
    int blocksDone = 0;
    for (query_idx = 0; query_idx < ref_n; query_idx++) {
      blocksDone += (ref_n + blockDim.x - 2 - query_idx) / blockDim.x;
      if (blocksDone > blockIdx.x) {
        break;
      }
    }
    // j (column) is given by multiplying the blocks needed for this i (row)
    // by the block size, plus offsets of i + 1 and the thread index
    int blocksPerQuery = (ref_n + blockDim.x - 2 - query_idx) / blockDim.x;
    ref_idx = query_idx + 1 + threadIdx.x +
              (blockIdx.x - (blocksDone - blocksPerQuery)) * blockDim.x;

    if (ref_idx < ref_n) {
      // Order of ref/query reversed here to give correct output order
      dist_idx = square_to_condensed(query_idx, ref_idx, ref_n);
    }
  } else {
    int blocksPerQuery = (ref_n + blockDim.x - 1) / blockDim.x;
    query_idx = blockIdx.x / blocksPerQuery;
    ref_idx = (blockIdx.x % blocksPerQuery) * blockDim.x + threadIdx.x;
    dist_idx = query_idx * ref_n + ref_idx;
  }
  __syncwarp();

  const uint64_t *ref_start = ref + ref_idx * ref_strides.sample_stride;
  const uint64_t *query_start = query + query_idx * query_strides.sample_stride;
  const float tolerance =
      __fdividef(5.0f, __int2float_rz(64 * ref_strides.sketchsize64));

  // Calculate Jaccard distances over k-mer lengths
  int kmer_used = 0;
  float xsum = 0.0f;
  float ysum = 0.0f;
  float xysum = 0.0f;
  float xsquaresum = 0.0f;
  float ysquaresum = 0.0f;
  bool stop = false;
  for (int kmer_idx = 0; kmer_idx < kmer_n; kmer_idx++) {
    // Copy query sketch into __shared__ mem
    // Uses all threads *in a single warp* to do the copy
    // NB there is no disadvantage vs using multiple warps, as they would have
    // to wait (see
    // https://stackoverflow.com/questions/15468059/copy-to-the-shared-memory-in-cuda)
    // NB for query these reads will be coalesced, but for ref they won't, as
    // can't coalesce both here (bin inner stride) and in jaccard (sample inner
    // stride)
    const uint64_t *query_ptr;
    extern __shared__ uint64_t query_shared[];
    int query_bin_strides;
    auto block = cooperative_groups::this_thread_block();
    __shared__ cuda::barrier<cuda::thread_scope::thread_scope_block> barrier;
    if (block.thread_rank() == 0) {
      init(&barrier, block.size()); // Friend function initializes barrier
    }
    block.sync();
    if (use_shared) {
      size_t sketch_bins = query_strides.bbits * query_strides.sketchsize64;
      size_t sketch_stride = query_strides.bin_stride;
      if (threadIdx.x < warp_size) {
        for (int lidx = threadIdx.x; lidx < sketch_bins; lidx += warp_size) {
          cuda::memcpy_async(query_shared + lidx,
                             query_start + (lidx * sketch_stride),
                             sizeof(uint64_t),
                             barrier);
        }
      }
      query_ptr = query_shared;
      query_bin_strides = 1;
    } else {
      query_ptr = query_start;
      query_bin_strides = query_strides.bin_stride;
    }
    barrier.arrive_and_wait();

    // Some threads at the end of the last block will have nothing to do
    // Need to have conditional here to avoid block on __syncthreads() above
    if (ref_idx < ref_n) {
      // Calculate Jaccard distance at current k-mer length
      float jaccard_obs = jaccard_dist(
          ref_start, query_ptr, ref_strides.sketchsize64, ref_strides.bbits,
          ref_strides.bin_stride, query_bin_strides);

      // Adjust for random matches
      float jaccard_expected =
          random_table[kmer_idx * random_strides.kmer_stride +
                       ref_idx_lookup[ref_idx] *
                           random_strides.cluster_inner_stride +
                       query_idx_lookup[query_idx] *
                           random_strides.cluster_outer_stride];
      float jaccard = observed_excess(jaccard_obs, jaccard_expected, 1.0f);
      // Stop regression if distances =~ 0
      if (jaccard < tolerance) {
        // Would normally break here, but gives no advantage on a GPU as causes
        // warp to diverge
        // As the thread blocks are used to load the query in, adding a break
        // would actually cause a stall. So just stop adding
        stop = true;
      } else if (!stop) {
        float y = __logf(jaccard);
        // printf("i:%d j:%d k:%d r:%f jac:%f y:%f\n", ref_idx, query_idx,
        // kmer_idx, jaccard_expected, jaccard_obs, y);

        // Running totals for regression
        kmer_used++;
        int kmer = kmers[kmer_idx];
        xsum += kmer;
        ysum += y;
        xysum += kmer * y;
        xsquaresum += kmer * kmer;
        ysquaresum += y * y;
      }
    }

    // Move to next k-mer length
    ref_start += ref_strides.kmer_stride;
    query_start += query_strides.kmer_stride;
  }

  if (ref_idx < ref_n) {
    // Run the regression, and store results in dists
    float fitted_dists[2];
    // Set diagonal if you wish to ignore diagonals (set them to max)
    if (max_diagonal && ref_idx == query_idx) {
      fitted_dists[0] = INFINITY;
      fitted_dists[1] = INFINITY;
    } else {
      simple_linear_regression(fitted_dists, xsum, ysum, xysum, xsquaresum,
                               ysquaresum, kmer_used, dist_col);
    }
    if (dist_col < 0) {
      dists[dist_idx] = fitted_dists[0];
      dists[dist_idx + dist_n] = fitted_dists[1];
    } else {
      dists[dist_idx] = fitted_dists[dist_col];
    }

    update_progress(dist_idx, dist_n, progress);
  }
}

/***************
 *			       *
 *	Host code  *
 *			       *
 ***************/

// Get the blockSize and blockCount for CUDA call
std::tuple<size_t, size_t> getBlockSize(const size_t ref_samples,
                                        const size_t query_samples,
                                        const size_t dist_rows,
                                        const bool self) {
  // Each block processes a single query. As max size is 512 threads
  // per block, may need multiple blocks (non-exact multiples lead
  // to some wasted computation in threads)
  // We take the next multiple of 32 that is larger than the number of
  // reference sketches, up to a maximum of 512
  size_t blockSize =
      std::min(256, 32 * static_cast<int>((ref_samples + 32 - 1) / 32));
  size_t blockCount = 0;
  if (self) {
    for (int i = 0; i < ref_samples; i++) {
      blockCount += (ref_samples + blockSize - 2 - i) / blockSize;
    }
  } else {
    size_t blocksPerQuery = (ref_samples + blockSize - 1) / blockSize;
    blockCount = blocksPerQuery * query_samples;
  }
  return (std::make_tuple(blockSize, blockCount));
}

// Writes a progress meter using the device int which keeps
// track of completed jobs
void reportDistProgress(progress_atomics& progress, long long dist_rows) {
  int now_completed = 0;
  float kern_progress = 0;
  if (dist_rows > progress_blocks) {
    while (now_completed < progress_blocks - 1) {
      if (PyErr_CheckSignals() != 0) {
        progress.set_kill();
        throw py::error_already_set();
      }
      int complete = progress.complete();
      if (complete > now_completed) {
        now_completed = complete;
        kern_progress = now_completed / (float)progress_blocks;
        fprintf(stderr, "%cProgress (GPU): %.1lf%%", 13, kern_progress * 100);
      } else {
        usleep(1000);
      }
    }
  }
}

// Initialise device and return info on its memory
std::tuple<size_t, size_t, size_t> initialise_device(const int device_id) {
  CUDA_CALL(cudaSetDevice(device_id));

  size_t mem_free = 0;
  size_t mem_total = 0;
  CUDA_CALL(cudaMemGetInfo(&mem_free, &mem_total));
  int shared_size = 0;
  CUDA_CALL(cudaDeviceGetAttribute(
      &shared_size, cudaDevAttrMaxSharedMemoryPerBlock, device_id));
  return (
      std::make_tuple(mem_free, mem_total, static_cast<size_t>(shared_size)));
}

std::tuple<bool, size_t> check_shared_size(const SketchStrides& strides, const size_t shared_size) {
  size_t sketch_size_bytes =
    strides.sketchsize64 * strides.bbits * sizeof(uint64_t);
  bool use_shared = true;
  if (sketch_size_bytes > shared_size) {
    std::cerr << "You are using a large sketch size, which may slow down "
                "computation on this device"
              << std::endl;
    std::cerr << "Reduce sketch size to "
              << std::floor(64 * shared_size /
                            (strides.bbits * sizeof(uint64_t)))
              << " or less for better performance" << std::endl;
    sketch_size_bytes = 0;
    use_shared = false;
  }
  return std::make_tuple(use_shared, sketch_size_bytes);
}

// Main function to run the distance calculations, reading/writing into
// device_arrays Cache preferences: Upper dist memory access is hard to predict,
// so try and cache as much as possible Query uses on-chip cache (__shared__) to
// store query sketch
std::vector<float> dispatchDists(std::vector<Reference> &ref_sketches,
                                 std::vector<Reference> &query_sketches,
                                 SketchStrides &ref_strides,
                                 SketchStrides &query_strides,
                                 const FlatRandom &flat_random,
                                 const std::vector<uint16_t> &ref_random_idx,
                                 const std::vector<uint16_t> &query_random_idx,
                                 const SketchSlice &sketch_subsample,
                                 const std::vector<size_t> &kmer_lengths,
                                 const bool self, const int cpu_threads,
                                 const size_t shared_size) {
  CUDA_CALL(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
  CUDA_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));

  // Progress meter
  progress_atomics progress;

  RandomStrides random_strides = std::get<0>(flat_random);
  long long dist_rows;
  if (self) {
    dist_rows =
        (sketch_subsample.ref_size * (sketch_subsample.ref_size - 1)) >> 1;
  } else {
    dist_rows = sketch_subsample.ref_size * sketch_subsample.query_size;
  }

  // Load memory onto device
  DeviceMemory device_arrays(ref_strides, query_strides, ref_sketches,
                             query_sketches, sketch_subsample, flat_random,
                             ref_random_idx, query_random_idx, kmer_lengths,
                             dist_rows, self, cpu_threads);

  bool use_shared;
  size_t shared_size_bytes;
  std::tie(use_shared, shared_size_bytes) =
    check_shared_size(query_strides, shared_size);

  size_t blockSize, blockCount;
  bool max_diagonal = false;
  int dist_col = -1;
  if (self) {
    std::tie(blockSize, blockCount) = getBlockSize(
        sketch_subsample.ref_size, sketch_subsample.ref_size, dist_rows, self);

    // Third argument is the size of __shared__ memory needed by a thread block
    // This is equal to the query sketch size in bytes (at a single k-mer
    // length)
    calculate_dists<<<blockCount, blockSize, shared_size_bytes>>>(
        self, device_arrays.ref_sketches(), sketch_subsample.ref_size,
        device_arrays.ref_sketches(), sketch_subsample.ref_size,
        device_arrays.kmers(), kmer_lengths.size(), device_arrays.dist_mat(),
        dist_rows, device_arrays.random_table(), device_arrays.ref_random(),
        device_arrays.ref_random(), ref_strides, ref_strides, random_strides,
        progress.get_ptrs(), use_shared, dist_col, max_diagonal);
  } else {
    std::tie(blockSize, blockCount) =
        getBlockSize(sketch_subsample.ref_size, sketch_subsample.query_size,
                     dist_rows, self);

    // Third argument is the size of __shared__ memory needed by a thread block
    // This is equal to the query sketch size in bytes (at a single k-mer
    // length)
    calculate_dists<<<blockCount, blockSize, shared_size_bytes>>>(
        self, device_arrays.ref_sketches(), sketch_subsample.ref_size,
        device_arrays.query_sketches(), sketch_subsample.query_size,
        device_arrays.kmers(), kmer_lengths.size(), device_arrays.dist_mat(),
        dist_rows, device_arrays.random_table(), device_arrays.ref_random(),
        device_arrays.query_random(), ref_strides, query_strides,
        random_strides, progress.get_ptrs(), use_shared, dist_col, max_diagonal);
  }

  // Check for error in kernel launch
  CUDA_CALL(cudaGetLastError());
  reportDistProgress(progress, dist_rows);
  fprintf(stderr, "%cProgress (GPU): 100.0%%\n", 13);

  // Copy results back to host
  CUDA_CALL(cudaDeviceSynchronize());
  std::vector<float> dist_results = device_arrays.read_dists();

  return (dist_results);
}


// Function which sparsifies distances on the fly. Distances are calculated in
// blocks, sorted and top k stored.
// NB cuda graph not needed as API calls faster than ops here
sparse_coo sparseDists(const dist_params params,
  const std::vector<std::vector<uint64_t>> &ref_sketches,
  const std::vector<SketchStrides> &ref_strides,
  const FlatRandom &flat_random,
  const std::vector<uint16_t> &ref_random_idx,
  const std::vector<size_t> &kmer_lengths,
  const int kNN,
  const size_t dist_col,
  const size_t samples_per_chunk,
  const size_t num_big_chunks,
  const int cpu_threads) {
  /*
   *
   *   Device setup
   *
   */
  CUDA_CALL(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
  CUDA_CALL(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));

  // Progress meter (not used for printing in this function)
  progress_atomics progress;

  // const parameters for functions
  const size_t n_chunks = ref_sketches.size();
  const bool self = false;
  const double total_blocks = static_cast<double>(n_chunks) * n_chunks;
  const size_t idx_blockSize = 64;
  const size_t copy_blockSize = 64;
  const int begin_sort_bit = 0;
  const int end_sort_bit = 8 * sizeof(float);

  bool use_shared;
  size_t shared_size_bytes;
  std::tie(use_shared, shared_size_bytes) =
    check_shared_size(ref_strides[0], params.shared_size);

  /*
   *
   *   Memory allocations
   *   All use the big_chunk size, for smaller chunks the end of the alloc
   *   is not used
   *
   */

  // Storage for results on host
  std::vector<float> host_dists(params.n_samples * kNN);
  std::vector<long> i_vec(params.n_samples * kNN);
  std::vector<long> j_vec(params.n_samples * kNN);

  //   sketch block 1
  //   sketch block 2
  //   sketch block 3
  //   sketch block 4
  device_array<uint64_t> block1(ref_sketches[0].size());
  device_array<uint64_t> block2(block1.size());
  device_array<uint64_t> block3(block1.size());
  device_array<uint64_t> block4(block1.size());
  //   dists (chunk size^2)
  //   dists idx (dists.size())
  device_array<float> dists((samples_per_chunk + 1) * (samples_per_chunk + 1));
  device_array<long> dists_idx(dists.size());

  // inner loop sorting
  //   sorted dists (dists.size())
  //   sorted dists idx (dists.size())
  //   offsets
  //   tmp space
  device_array<float> sorted_dists(dists.size());
  device_array<long> sorted_dists_idx(dists.size());
  device_array<int> dist_partitions((samples_per_chunk + 1) + 1);
  device_array<void> dist_sort_tmp;

  // inner loop kNN pick (copy_top_k)
  //   sorted dists
  //   sorted dists idx
  device_array<float> all_sorted_dists(kNN * n_chunks * (samples_per_chunk + 1));
  device_array<long> all_sorted_dists_idx(all_sorted_dists.size());

  // outer loop sorting
  //   doubly sorted dists (staging for final results, copied to host)
  //   doubly sorted dists idx (as above)
  //   offsets
  //   tmp space
  device_array<float> doubly_sorted_dists(all_sorted_dists.size());
  device_array<long> doubly_sorted_dists_idx(all_sorted_dists.size());
  std::vector<int> host_partitions;
  for (int partition_idx = 0; partition_idx < dist_partitions.size(); ++partition_idx) {
    host_partitions.push_back(partition_idx * kNN * n_chunks);
  }
  device_array<int> second_sort_partitions(host_partitions);
  host_partitions.clear();
  device_array<void> outer_dist_sort_tmp;

  // out loop kNN pick (copy_top_k)
  //   final dists (staging for final results, copied to host)
  //   final dists idx (as above)
  device_array<float> final_dists(kNN * (samples_per_chunk + 1));
  device_array<long> final_dists_idx(kNN * (samples_per_chunk + 1));

  //   kmers
  const std::vector<int> kmer_ints(kmer_lengths.begin(), kmer_lengths.end());
  device_array<int> kmers(kmer_ints);

  //   random
  RandomStrides random_strides = std::get<0>(flat_random);
  device_array<float> random_table(std::get<1>(flat_random));
  device_array<uint16_t> random_idx(ref_random_idx);

  /*
   *
   *   Memory setup
   *   Set i_vec on host
   *   Register sketches in host memory
   *   Copy first sketch blocks to device
   *
   */

  // Set i_vec (it just needs to be 0, 0, 0..., 1, 1, 1...)
  #pragma omp parallel for num_threads(cpu_threads)
  for (size_t sample_idx = 0; sample_idx < params.n_samples; ++sample_idx) {
    std::fill_n(i_vec.begin() + sample_idx * kNN, kNN, sample_idx);
  }

  // Register sketches on host so they can be copied async
  for (size_t chunk_idx = 0; chunk_idx < n_chunks; ++chunk_idx) {
    CUDA_CALL(cudaHostRegister((void *)ref_sketches[chunk_idx].data(),
                               ref_sketches[chunk_idx].size() * sizeof(uint64_t),
                               cudaHostRegisterReadOnly));
  }

  // Four CUDA streams used in loops
  cuda_stream dist_stream, idx_stream, mem_stream, sort_stream;

  // Copy first set of refs in. Block 4 is used to store this permanently to stop
  // two copies being needed when entering a new row
  block1.set_array_async(ref_sketches[0].data(), ref_sketches[0].size(), mem_stream.stream());
  block4.set_array_async(block1.data(), block1.size(), mem_stream.stream());

  /*
   *
   *   Loop over chunks
   *   Outer loop over chunks of rows
   *   Inner loop over chunks of columns
   *
   */
  fprintf(stderr, "%cProgress (GPU): %.1lf%%", 13, 0.0f);

  // OUTER LOOP over n_chunks lots of refs
  size_t row_offset = 0;
  for (size_t row_chunk_idx = 0; row_chunk_idx < n_chunks; ++row_chunk_idx) {
    size_t row_samples = samples_per_chunk + (row_chunk_idx < num_big_chunks ? 1 : 0);
    size_t col_offset = 0;

    //  INNER LOOP over n_chunks lots of queries
    for (size_t col_chunk_idx = 0; col_chunk_idx < n_chunks; ++col_chunk_idx) {
      size_t col_samples = samples_per_chunk + (col_chunk_idx < num_big_chunks ? 1 : 0);
      size_t dist_rows = row_samples * col_samples;

      // Check for interrupts
      if (PyErr_CheckSignals() != 0) {
        throw py::error_already_set();
      }

      // Only need to set new sort partitions if moving to a smaller chunk
      if (host_partitions.size() != col_samples + 1) {
        host_partitions.clear();
        for (int partition_idx = 0; partition_idx < col_samples + 1; ++partition_idx) {
          host_partitions.push_back(partition_idx * col_samples);
        }
        dist_partitions.set_array_async(host_partitions.data(), host_partitions.size(), sort_stream.stream());
      }

      //    (stream 1 async) Run dists on 1 vs 2
      size_t blockSize, blockCount;
      std::tie(blockSize, blockCount) =
        getBlockSize(row_samples, col_samples, dist_rows, self);
      uint64_t* query_ptr = col_chunk_idx == 0 ? block4.data() : block2.data();
      bool max_diagonal = col_chunk_idx == row_chunk_idx;
      mem_stream.sync();
      // NB in calculate_dists ref idx changes fastest (so should be the column)
      // so ref and query are 'backwards'
      calculate_dists<<<blockCount, blockSize, shared_size_bytes, dist_stream.stream()>>>(
        self,
        query_ptr, col_samples,
        block1.data(), row_samples,
        kmers.data(), kmers.size(),
        dists.data(), dist_rows,
        random_table.data(), random_idx.data() + col_offset, random_idx.data() + row_offset,
        ref_strides[col_chunk_idx], ref_strides[row_chunk_idx],
        random_strides, progress.get_ptrs(), use_shared, dist_col, max_diagonal
      );

      //    (stream 2 async) Load next into block 3
      //    swap ptrs for block 2 <-> 3
      if (col_chunk_idx + 1 < n_chunks) {
        block3.set_array_async(ref_sketches[col_chunk_idx + 1].data(), ref_sketches[col_chunk_idx + 1].size(), mem_stream.stream());
        block2.swap(block3);
      } else if (row_chunk_idx + 1 < n_chunks) {
        block3.set_array_async(ref_sketches[row_chunk_idx + 1].data(), ref_sketches[row_chunk_idx + 1].size(), mem_stream.stream());
        block1.swap(block3);
      }

      //    (stream 3 async) Set dist idx via kernel
      const size_t idx_blockCount = (dist_rows + idx_blockSize - 1) / idx_blockSize;
      set_idx<<<idx_blockCount, idx_blockSize, 0, idx_stream.stream()>>>(
        dists_idx.data(), row_samples, col_samples, col_offset
      );

      //    (stream 4) cub::DeviceSegmentedRadixSort::SortPairs on dists, dists idx -> sorted dists, sorted dists idx
      const int num_items = dist_rows;
      const int num_segments = host_partitions.size() - 1;
      int *d_offsets = dist_partitions.data();
      float *d_keys_in = dists.data();
      float *d_keys_out = sorted_dists.data();
      long *d_values_in = dists_idx.data();
      long *d_values_out = sorted_dists_idx.data();
      // Determine temporary device storage requirements (first run only)
      if (row_chunk_idx == 0 && col_chunk_idx == 0) {
        size_t temp_storage_bytes = 0;
        cub::DeviceSegmentedRadixSort::SortPairs(dist_sort_tmp.data(), temp_storage_bytes,
            d_keys_in, d_keys_out, d_values_in, d_values_out,
            num_items, num_segments, d_offsets, d_offsets + 1,
            begin_sort_bit, end_sort_bit, sort_stream.stream());
        dist_sort_tmp.set_size(temp_storage_bytes);
      }
      // Run sorting operation
      dist_stream.sync();
      idx_stream.sync();
      size_t temp_storage_bytes = dist_sort_tmp.size();
      cub::DeviceSegmentedRadixSort::SortPairs(dist_sort_tmp.data(), temp_storage_bytes,
          d_keys_in, d_keys_out, d_values_in, d_values_out,
          num_items, num_segments, d_offsets, d_offsets + 1,
          begin_sort_bit, end_sort_bit, sort_stream.stream());

      //    (stream 4) D->D copy of top kNN dists to start of dists
      const bool second_sort = false;
      const size_t dist_out_size = kNN * num_segments;
      const size_t copy_blockCount = (dist_out_size + copy_blockSize - 1) / copy_blockSize;
      copy_top_k<<<copy_blockCount, copy_blockSize, 0, sort_stream.stream()>>>(
        sorted_dists.data(), sorted_dists_idx.data(), all_sorted_dists.data(), all_sorted_dists_idx.data(),
        col_samples, dist_out_size, kNN, col_chunk_idx, n_chunks, second_sort
      );

      // Update progress
      col_offset += col_samples;
      const size_t blocks_done = row_chunk_idx * n_chunks + col_chunk_idx + 1;
      fprintf(stderr, "%cProgress (GPU): %.1lf%%", 13, 100 * blocks_done / total_blocks);
    }

    const int num_items = all_sorted_dists.size();
    const int num_segments = row_samples;
    int *d_offsets = second_sort_partitions.data();
    float *d_keys_in = all_sorted_dists.data();
    float *d_keys_out = doubly_sorted_dists.data();
    long *d_values_in = all_sorted_dists_idx.data();
    long *d_values_out = doubly_sorted_dists_idx.data();
    // Determine temporary device storage requirements
    if (row_chunk_idx == 0) {
      size_t temp_storage_bytes = 0;
      cub::DeviceSegmentedRadixSort::SortPairs(outer_dist_sort_tmp.data(), temp_storage_bytes,
          d_keys_in, d_keys_out, d_values_in, d_values_out,
          num_items, num_segments, d_offsets, d_offsets + 1,
          begin_sort_bit, end_sort_bit, sort_stream.stream());
      outer_dist_sort_tmp.set_size(temp_storage_bytes);
    }
    // Run sorting operation
    size_t temp_storage_bytes = outer_dist_sort_tmp.size();
    cub::DeviceSegmentedRadixSort::SortPairs(outer_dist_sort_tmp.data(), temp_storage_bytes,
        d_keys_in, d_keys_out, d_values_in, d_values_out,
        num_items, num_segments, d_offsets, d_offsets + 1,
        begin_sort_bit, end_sort_bit, sort_stream.stream());

    // take top kNN
    const bool second_sort = true;
    const int block_offset = 0;
    const size_t dist_out_size = kNN * num_segments;
    const size_t copy_blockCount = (dist_out_size + copy_blockSize - 1) / copy_blockSize;
    copy_top_k<<<copy_blockCount, copy_blockSize, 0, sort_stream.stream()>>>(
      doubly_sorted_dists.data(), doubly_sorted_dists_idx.data(), final_dists.data(), final_dists_idx.data(),
      kNN * n_chunks, dist_out_size, kNN, block_offset, n_chunks, second_sort
    );
    // Copy chunk of results back to host
    final_dists.get_array_async(host_dists.data() + row_offset * kNN, dist_out_size, sort_stream.stream());
    final_dists_idx.get_array_async(j_vec.data() + row_offset * kNN, dist_out_size, sort_stream.stream());

    row_offset += row_samples;
  }
  fprintf(stderr, "%cProgress (GPU): 100.0%%\n", 13);
  CUDA_CALL(cudaDeviceSynchronize());

  // Unregister host memory
  for (size_t chunk_idx = 0; chunk_idx < n_chunks; ++chunk_idx) {
    CUDA_CALL(cudaHostUnregister((void *)ref_sketches[chunk_idx].data()));
  }

  return (std::make_tuple(i_vec, j_vec, host_dists));
}
