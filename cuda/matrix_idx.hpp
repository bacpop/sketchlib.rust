/*
 *
 * matrix_idx.hpp
 * long/square index conversion
 *
 */
#pragma once

#include <cstddef> // size_t
#include <cmath>   // floor/sqrt
#include <cassert> // assert

template <class T>
inline size_t rows_to_samples(const T &longMat)
{
  return 0.5 * (1 + sqrt(1 + 8 * (longMat.rows())));
}

#ifdef __NVCC__
__host__ __device__
#endif
    inline long
    calc_row_idx(const long long k, const long n)
{
#ifndef __CUDA_ARCH__
  return n - 2 - std::floor(std::sqrt(static_cast<double>(-8 * k + 4 * n * (n - 1) - 7)) / 2 - 0.5);
#else
  // __ll2float_rn() casts long long to float, rounding to nearest
  return n - 2 - floor(__dsqrt_rn(__ll2double_rz(-8 * k + 4 * n * (n - 1) - 7)) / 2 - 0.5);
#endif
}

#ifdef __NVCC__
__host__ __device__
#endif
    inline long
    calc_col_idx(const long long k, const long i, const long n)
{
  return k + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2;
}

#ifdef __NVCC__
__host__ __device__
#endif
    inline long long
    square_to_condensed(long i, long j, long n)
{
  assert(j > i);
  return (n * i - ((i * (i + 1)) >> 1) + j - 1 - i);
}
