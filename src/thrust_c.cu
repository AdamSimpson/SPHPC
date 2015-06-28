#define restrict __restrict

#include "thrust_c.h"
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <iostream>

struct LessThan {
  LessThan(double min): x_min(min) {}

  double x_min;

  __host__ __device__
  bool operator()(const double x) {
    return (x < x_min);
  }
};

struct LessThanOrEqual {
  LessThanOrEqual(double min): x_min(min) {}

  double x_min;

  __host__ __device__
  bool operator()(const double x) {
    return (x <= x_min);
  }
};

struct GreaterThan {
  GreaterThan(double max): x_max(max) {}

  double x_max;

  __host__ __device__
  bool operator()(const double x) {
    return (x > x_max);
  }
};

struct GreaterThanOrEqual {
  GreaterThanOrEqual(double max): x_max(max) {}

  double x_max;

  __host__ __device__
  bool operator()(const double x) {
    return (x > x_max);
  }
};

struct OutsideBounds {
  OutsideBounds(double min, double max): x_min(min), x_max(max) {}

  double x_min, x_max;

  __host__ __device__
  bool operator()(const double x) {
    return (x < x_min || x > x_max);
  }
};

extern "C" void SortByKey(unsigned int *const keys,
                          unsigned int *const values,
                          const int count) {
  thrust::sort_by_key(thrust::device, keys, keys + count, values);

}

extern "C" void FindLowerBounds(const unsigned int *const hash_values,
                                const int value_count,
                                const int search_count,
                                unsigned int *const start_indices) {

  // Why do I need this when i'm using device execution policy?
  thrust::device_ptr<const unsigned int> d_hash_values = thrust::device_pointer_cast(hash_values);
  thrust::device_ptr<unsigned int> d_start_indices = thrust::device_pointer_cast(start_indices);

  thrust::lower_bound(thrust::device,
                      d_hash_values,
                      d_hash_values + value_count,
                      thrust::make_counting_iterator((unsigned int)0),
                      thrust::make_counting_iterator((unsigned int)search_count),
                      d_start_indices);
}

extern "C" void FindUpperBounds(const unsigned int *const hash_values,
                                const int value_count,
                                const int search_count,
                                unsigned int *const end_indices) {

  // Why do I need this when i'm using device execution policy?
  thrust::device_ptr<const unsigned int> d_hash_values = thrust::device_pointer_cast(hash_values);
  thrust::device_ptr<unsigned int> d_end_indices = thrust::device_pointer_cast(end_indices);

  thrust::upper_bound(thrust::device,
                      hash_values,
                      hash_values + value_count,
                      thrust::make_counting_iterator((unsigned int)0),
                      thrust::make_counting_iterator((unsigned int)search_count),
                      d_end_indices);
}

// C wrapper function for thrust copy_if with LessThan predicate
extern "C" void CopyIfLessThan(const double min,
                               const int *const input,
                               const int input_count,
                               const double *const stencil,
                               int *const output,
                               int *const output_count) {

  int *end_pointer = thrust::copy_if(thrust::device,
                                     input,
                                     input + input_count,
                                     stencil,
                                     output,
                                     LessThan(min));
  *output_count = end_pointer - output;
}

// C wrapper function for thrust copy_if with LessThanOrEqual predicate
extern "C" void CopyIfLessThanOrEqual(const double min,
                                      const int *const input,
                                      const int input_count,
                                      const double *const stencil,
                                      int *const output,
                                      int *const output_count) {

  int *end_pointer = thrust::copy_if(thrust::device,
                                     input,
                                     input + input_count,
                                     stencil,
                                     output,
                                     LessThanOrEqual(min));
  *output_count = end_pointer - output;
}

// C wrapper function for thrust copy_if with GreaterThan predicate
extern "C" void CopyIfGreaterThan(const double max,
                                  const int *const input,
                                  const int input_count,
                                  const double *const stencil,
                                  int *const output,
                                  int *const output_count) {

  int *end_pointer = thrust::copy_if(thrust::device,
                                     input,
                                     input + input_count,
                                     stencil,
                                     output,
                                     GreaterThan(max));
  *output_count = end_pointer - output;
}

// C wrapper function for thrust copy_if with GreaterThan predicate
extern "C" void CopyIfGreaterThanOrEqual(const double max,
                                         const int *const input,
                                         const int input_count,
                                         const double *const stencil,
                                         int *const output,
                                         int *const output_count) {

  int *end_pointer = thrust::copy_if(thrust::device,
                                     input,
                                     input + input_count,
                                     stencil,
                                     output,
                                     GreaterThanOrEqual(max));
  *output_count = end_pointer - output;
}

// C wrapper function for thrust remove_if with OutsideBounds predicate
extern "C" void RemoveIfOutsideBounds(const double min, const double max,
                                      int *const input,
                                      const int input_count,
                                      const double *const stencil) {

  thrust::remove_if(thrust::device,
                    input,
                    input + input_count,
                    stencil,
                    OutsideBounds(min, max));
}
