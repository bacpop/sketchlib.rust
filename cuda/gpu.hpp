/*
 *
 * gpu.hpp
 * functions using CUDA
 *
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "device_reads.cuh"
#include "dist/matrix_types.hpp"

static const int warp_size = 32;

// Small struct used in cuda_dists_init
struct ALIGN(16) dist_params {
  bool self;
  SketchStrides ref_strides;
  SketchStrides query_strides;
  long long dist_rows;
  long n_samples;
  size_t shared_size;
  size_t mem_free;
  size_t mem_total;
};

// defined in dist.cu
std::tuple<size_t, size_t, size_t> initialise_device(const int device_id);

std::vector<uint64_t> flatten_by_bins(const std::vector<Reference> &sketches,
                                      const std::vector<size_t> &kmer_lengths,
                                      SketchStrides &strides,
                                      const size_t start_sample_idx,
                                      const size_t end_sample_idx,
                                      const int cpu_threads = 1);

std::vector<uint64_t>
flatten_by_samples(const std::vector<Reference> &sketches,
                   const std::vector<size_t> &kmer_lengths,
                   SketchStrides &strides, const size_t start_sample_idx,
                   const size_t end_sample_idx, const int cpu_threads = 1);

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
                                 const size_t shared_size);

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
  const int cpu_threads);


