#pragma once

#if !defined(NO_CUDA)
#include "cuda_runtime.h"
#endif

#include <cstddef>

/**
   Simple array like class to handle pool of managed memory
 **/
namespace sim {

template <typename T>
class Array {
public:
  Array(std::size_t capacity) : capacity_{capacity}, size_{0} {
    #if defined(NO_CUDA)
    data_ = new T [capacity];
    #else
    cudaMallocManaged(&data_, capacity_*sizeof(T));
    #endif
  }

  ~Array() {
    #if defined(NO_CUDA)
    free(data_);
    #else
    cudaFree(data_);
    #endif
  }

  /**
     Copying/moving the array is bad as it can free the data prematurely
     Copy the raw pointer instead using data() for tmp variables
   **/
  Array(const Array&)            =delete;
  Array& operator=(const Array&) =delete;
  Array(Array&&) noexcept        =delete;
  Array& operator=(Array&&)      =delete;

  std::size_t size() const {
    return size_;
  }

  std::size_t capacity() const {
    return capacity_;
  }

  std::size_t available() const {
    return capacity() - size();
  }

  void push_back(const T& value) {
    data_[size_++] = value;
  }

  void push_back(const T* values_ptr, const size_t push_count) {
    for(std::size_t i=0; i<push_count; i++)
      this->push_back(values_ptr[i]);
  }

  void pop_back() {
    size_--;
  }

  void pop_back(std::size_t pop_count) {
    size_ -= pop_count;
  }

  T* data() {
    return this->data_;
  }

  T* data() const {
    return this->data_;
  }

  T& operator[] (const std::size_t index) {
    return data_[index];
  }

  const T& operator[] (const std::size_t index) const {
    return data_[index];
  }

private:
  T* data_;
  std::size_t capacity_;
  std::size_t size_;
};
}
