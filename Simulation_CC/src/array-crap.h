#pragma once

#include "thrust/host_vector.h"
#include "thrust/device_vector.h"
#include "thrust/execution_policy.h"
#include <algorithm>

/**
   Array class to handle host/device thrust vectors
   A subet of thrust vector functions are implemented
 **/
namespace sim {

template <typename T>
class Vector {
public:
  Vector() {};
  Vector(std::size_t n) : d_vec{n}, h_vec{n} {};
  Vector(std::size_t n, const T &value) : d_vec{n, value}, h_vec{n, value} {};

  void resize(const std::size_t new_size) {
    d_vec.resize(new_size);
    h_vec.resize(new_size);
  }

  const std::size_t size() const {
    assert(h_vec.size() == d_vec.size());
    return h_vec.size();
  }

  const std::size_t max_size() const {
    std::size_t h_size = h_vec.size();
    std::size_t d_size = d_vec.size();

    return std::min(h_size, d_size);
  }

  void reserve(std::size_t new_capacity) {
    d_vec.reserve(new_capacity);
    h_vec.reserve(new_capacity);
  }

  thrust::device_vector<T>& host() {
    return h_vec;
  }

  thrust::host_vector<T>& device() {
    return d_vec;
  }

  const thrust::device_vector<T>& host() const {
    return h_vec;
  }

  const thrust::host_vector<T>& device() const {
    return d_vec;
  }

  T* host_data() {
    return h_vec.data();
  }

  T* device_data() {
    return d_vec.data();
  }

  T* data(const thrust::execution_policy<thrust::host>) {
    return this->host_data();
  }

  const T* host_data() const {
    return h_vec.data();
  }

  const T* device_data() const {
    return d_vec.data();
  }

  template <typename PolicyType>
  const T* data(const PolicyType policy) const {
    if(policy == thrust::host)
      return this->host_data();
    else if(policy == thrust::device)
      return this->device_data();
  }

  template <typename PolicyType>
  void push_back(const T& x, const PolicyType policy) {
    if(policy == thrust::host)
      h_vec.push_back(x);
    else if(policy == thrust::device)
      d_vec.push_back(x);
  }

private:
  thrust::device_vector<T> d_vec;
  thrust::host_vector<T> h_vec;
};
}
