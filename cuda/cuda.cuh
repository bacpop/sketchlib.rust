#pragma once

#include <stdexcept>
#include <stdio.h>

#include <type_traits>
static_assert(__CUDACC_VER_MAJOR__ >= 11, "CUDA >=11.0 required");

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#include "align.hpp"

const int progress_blocks = 1000; // Update 1000

static void HandleCUDAError(const char *file, int line,
                            cudaError_t status = cudaGetLastError()) {
#ifdef _DEBUG
  cudaDeviceSynchronize();
#endif

  if (status != cudaSuccess || (status = cudaGetLastError()) != cudaSuccess) {
    if (status == cudaErrorUnknown) {
      printf("%s(%i) An Unknown CUDA Error Occurred :(\n", file, line);
    }
    printf("%s(%i) CUDA Error Occurred;\n%s\n", file, line,
           cudaGetErrorString(status));
    throw std::runtime_error("CUDA error");
  }
}

#define CUDA_CALL(err) (HandleCUDAError(__FILE__, __LINE__, err))
#define CUDA_CALL_NOTHROW( err ) (err)

struct ALIGN(16) progress_ptrs {
  volatile int *blocks_complete;
  bool *kill_kernel;
};

class progress_atomics {
public:
  progress_atomics() {
    CUDA_CALL(cudaMallocManaged((void**)&(managed_ptrs.blocks_complete), sizeof(int)));
    CUDA_CALL(cudaMallocManaged((void**)&(managed_ptrs.kill_kernel), sizeof(bool)));
    *(managed_ptrs.blocks_complete) = 0;
    *(managed_ptrs.kill_kernel) = false;
  }

  ~progress_atomics() {
    CUDA_CALL(cudaFree((void *)managed_ptrs.blocks_complete));
    CUDA_CALL(cudaFree(managed_ptrs.kill_kernel));
  }

  void set_kill() {
    *(managed_ptrs.kill_kernel) = true;
  }

  progress_ptrs get_ptrs() {
    return managed_ptrs;
  }

  __host__
  int complete() {
    return *(managed_ptrs.blocks_complete);
  }

private:
  progress_atomics(const progress_atomics &) = delete;
  progress_atomics(progress_atomics &&) = delete;

  progress_ptrs managed_ptrs;
};

// Use atomic add to update a counter, so progress works regardless of
// dispatch order
__device__ inline void update_progress(long long dist_idx, long long dist_n,
                                       progress_ptrs &progress) {
  // Progress indicator
  if (dist_idx % (dist_n / progress_blocks) == 0) {
    if (*(progress.kill_kernel)) {
      __trap();
    }
    atomicAdd((int *)progress.blocks_complete, 1);
    __threadfence_system();
  }
}
