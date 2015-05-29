#define restrict __restrict

#include "thrust_c.h"
#include <thrust/execution_policy.h>
#include <thrust/sort.h>
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
  thrust::sort_by_key(thrust::host, keys, keys + count, values);

}

// C wrapper function for thrust copy_if with LessThan predicate
extern "C" void CopyIfLessThan(const double min,
                               const int *const input,
                               const int input_count,
                               const double *const stencil,
                               int *const output,
                               int *const output_count) {

  int *end_pointer = thrust::copy_if(thrust::host,
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

  int *end_pointer = thrust::copy_if(thrust::host,
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

  int *end_pointer = thrust::copy_if(thrust::host,
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

  int *end_pointer = thrust::copy_if(thrust::host,
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

  thrust::remove_if(thrust::host,
                    input,
                    input + input_count,
                    stencil,
                    OutsideBounds(min, max));
}
