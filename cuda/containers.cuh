#pragma once

#include <vector>

#include "cuda.cuh"

template <typename T> class device_array {
public:
  // Default constructor
  device_array() : data_(nullptr), size_(0) {}

  // Constructor to allocate empty memory
  device_array(const size_t size) : size_(size) {
    CUDA_CALL(cudaMalloc((void **)&data_, size_ * sizeof(T)));
    CUDA_CALL(cudaMemset(data_, 0, size_ * sizeof(T)));
  }

  // Constructor from vector
  device_array(const std::vector<T> &data) : size_(data.size()) {
    CUDA_CALL(cudaMalloc((void **)&data_, size_ * sizeof(T)));
    CUDA_CALL(
        cudaMemcpy(data_, data.data(), size_ * sizeof(T), cudaMemcpyDefault));
  }

  // Copy
  device_array(const device_array &other) : size_(other.size_) {
    CUDA_CALL(cudaMalloc((void **)&data_, size_ * sizeof(T)));
    CUDA_CALL(
        cudaMemcpy(data_, other.data_, size_ * sizeof(T), cudaMemcpyDefault));
  }

  // Copy assign
  device_array &operator=(const device_array &other) {
    if (this != &other) {
      size_ = other.size_;
      CUDA_CALL(cudaFree(data_));
      CUDA_CALL(cudaMalloc((void **)&data_, size_ * sizeof(T)));
      CUDA_CALL(
          cudaMemcpy(data_, other.data_, size_ * sizeof(T), cudaMemcpyDefault));
    }
    return *this;
  }

  // Move
  device_array(device_array &&other) : data_(nullptr), size_(0) {
    data_ = other.data_;
    size_ = other.size_;
    other.data_ = nullptr;
    other.size_ = 0;
  }

  // Move assign
  device_array &operator=(device_array &&other) {
    if (this != &other) {
      CUDA_CALL(cudaFree(data_));
      data_ = other.data_;
      size_ = other.size_;
      other.data_ = nullptr;
      other.size_ = 0;
    }
    return *this;
  }

  ~device_array() { CUDA_CALL_NOTHROW(cudaFree(data_)); }

  void swap(device_array& other) {
    T* data_tmp = other.data_;
    size_t size_tmp = other.size_;
    other.data_ = data_;
    other.size_ = size_;
    data_ = data_tmp;
    size_ = size_tmp;
  }

  void get_array(std::vector<T> &dst) const {
    CUDA_CALL(cudaMemcpy(dst.data(), data_, dst.size() * sizeof(T),
                         cudaMemcpyDefault));
  }

  void set_array(const std::vector<T> &src) {
    size_ = src.size();
    CUDA_CALL(
        cudaMemcpy(data_, src.data(), size_ * sizeof(T), cudaMemcpyDefault));
  }

  void set_array(const T *src) {
    CUDA_CALL(cudaMemcpy(data_, src, size_ * sizeof(T), cudaMemcpyDefault));
  }

  // Special functions for working with streams/graphs
  // (prefer above functions unless using these)
  void get_array_async(T *value, const size_t size, cudaStream_t stream) const {
    if (size > size_) {
      throw std::runtime_error("Size of host memory too big to copy");
    }
    CUDA_CALL(cudaMemcpyAsync((void *)value, data_, sizeof(T) * size,
                              cudaMemcpyDefault, stream));
  }

  void set_array_async(const T *value, const size_t size, cudaStream_t stream) {
    if (size > size_) {
      throw std::runtime_error("Size of host memory too big to copy");
    }
    CUDA_CALL(cudaMemcpyAsync(data_, (void *)value, sizeof(T) * size,
                              cudaMemcpyDefault, stream));
  }

  T *data() { return data_; }

  size_t size() const { return size_; }

private:
  T *data_;
  size_t size_;
};

// Specialisation of the above for void* memory needed by some cub functions
// Construct once and use set_size() to modify
template <> class device_array<void> {
public:
  // Default constructor
  device_array() : data_(nullptr), size_(0) {}
  // Constructor to allocate empty memory
  device_array(const size_t size) : size_(size) {
    if (size_ > 0) {
      CUDA_CALL(cudaMalloc((void **)&data_, size_));
    }
  }

  ~device_array() { CUDA_CALL_NOTHROW(cudaFree(data_)); }

  void set_size(size_t size) {
    if (size != size_) {
      size_ = size;
      CUDA_CALL(cudaFree(data_));
      if (size_ > 0) {
        CUDA_CALL(cudaMalloc((void **)&data_, size_));
      } else {
        data_ = nullptr;
      }
    }
  }

  void *data() { return data_; }

  size_t size() const { return size_; }

private:
  device_array(const device_array<void> &) = delete;
  device_array(device_array<void> &&) = delete;

  void *data_;
  size_t size_;
};

class cuda_stream {
public:
  cuda_stream() { CUDA_CALL(cudaStreamCreate(&stream_)); }

  ~cuda_stream() {
    if (stream_ != nullptr) {
      CUDA_CALL_NOTHROW(cudaStreamDestroy(stream_));
    }
  }

  cudaStream_t stream() { return stream_; }

  void sync() { CUDA_CALL(cudaStreamSynchronize(stream_)); }

private:
  // Delete copy and move
  cuda_stream(const cuda_stream &) = delete;
  cuda_stream(cuda_stream &&) = delete;

  cudaStream_t stream_;
};
