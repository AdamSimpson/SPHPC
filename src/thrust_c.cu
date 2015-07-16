#define restrict __restrict

#include "thrust_c.h"
#include <thrust/execution_policy.h>
#include <thrust/version.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/remove.h>
#include <thrust/binary_search.h>
#include <thrust/device_ptr.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <iostream>
#include <thrust/system_error.h>
//#include "cuda_runtime.h"

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

struct RemoveOutsideBounds {
  RemoveOutsideBounds(double min, double max): x_min(min), x_max(max) {}

  double x_min, x_max;

  template <typename T>
  __host__ __device__
  bool operator()(T &tuple) const {
    const double x_star = thrust::get<0>(tuple);
    return (x_star < x_min || x_star > x_max);
  }
};

extern "C" void SortByKey(unsigned int *const keys,
                          unsigned int *const values,
                          const int count,
                          const void *const cuda_stream) {
  thrust::sort_by_key(thrust::cuda::par.on((cudaStream_t)cuda_stream), keys, keys + count, values);
}

extern "C" void FindLowerBounds(const unsigned int *const hash_values,
                                const int value_count,
                                const int search_count,
                                unsigned int *const start_indices,
                                const void *const cuda_stream) {

  thrust::lower_bound(thrust::cuda::par.on((cudaStream_t)cuda_stream),
                      hash_values,
                      hash_values + value_count,
                      thrust::make_counting_iterator((unsigned int)0),
                      thrust::make_counting_iterator((unsigned int)search_count),
                      start_indices);
}

extern "C" void FindUpperBounds(const unsigned int *const hash_values,
                                const int value_count,
                                const int search_count,
                                unsigned int *const end_indices,
                                const void *const cuda_stream) {

  thrust::upper_bound(thrust::cuda::par.on((cudaStream_t)cuda_stream),
                      hash_values,
                      hash_values + value_count,
                      thrust::make_counting_iterator((unsigned int)0),
                      thrust::make_counting_iterator((unsigned int)search_count),
                      end_indices);

}

// C wrapper function for thrust copy_if with LessThan predicate
extern "C" void CopyIfLessThan(const double min,
                               const int *const input,
                               const int input_count,
                               const double *const stencil,
                               int *const output,
                               int *const output_count,
                               const void *const cuda_stream) {

  int *end_pointer = thrust::copy_if(thrust::cuda::par.on((cudaStream_t)cuda_stream),
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
                                      int *const output_count,
                                      const void *const cuda_stream) {

  int *end_pointer = thrust::copy_if(thrust::cuda::par.on((cudaStream_t)cuda_stream),
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
                                  int *const output_count,
                                  const void *const cuda_stream) {

  int *end_pointer = thrust::copy_if(thrust::cuda::par.on((cudaStream_t)cuda_stream),
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
                                         int *const output_count,
                                         const void *const cuda_stream) {

  int *end_pointer = thrust::copy_if(thrust::cuda::par.on((cudaStream_t)cuda_stream),
                                     input,
                                     input + input_count,
                                     stencil,
                                     output,
                                     GreaterThanOrEqual(max));
  *output_count = end_pointer - output;
}

extern "C" void SortParticlesByKey(unsigned int *const keys,
                                 const int count,
                                 double *x_star, double *y_star, double *z_star,
                                 double *x, double *y, double *z,
                                 double *v_x, double *v_y, double *v_z,
                                 const void *const cuda_stream) {

  thrust::sort_by_key(thrust::cuda::par.on((cudaStream_t)cuda_stream),
                 keys, keys+count,
                 thrust::make_zip_iterator(thrust::make_tuple(x_star, y_star, z_star,
                                                              x, y, z, v_x, v_y, v_z))
                );
}

// C wrapper function for thrust remove_if with OutsideBounds predicate
extern "C" void RemoveIfOutsideBounds(const double min, const double max,
                                      int *const input,
                                      const int input_count,
                                      const double *const stencil,
                                      const void *const cuda_stream) {

  int *end_pointer = thrust::remove_if(thrust::cuda::par.on((cudaStream_t)cuda_stream),
                       input,
                       input + input_count,
                       stencil,
                       OutsideBounds(min, max));
}
// Rename me
extern "C" void RemoveIfOutsideBounds2(const double min, const double max,
                                      int *const input,
                                      const int input_count,
                                      const void *const cuda_stream) {

  int *end_pointer = thrust::remove_if(thrust::cuda::par.on((cudaStream_t)cuda_stream),
                       input,
                       input + input_count,
                       OutsideBounds(min, max));
}

// Rename me
extern "C" void RemoveIfOutsideBounds3(const double min, const double max,
                                      const int count,
                                      double *x_star, double *y_star, double *z_star,
                                      double *x, double *y, double *z,
                                      double *v_x, double *v_y, double *v_z,
                                      const void *const cuda_stream) {
/*  typedef thrust::tuple< double*, double*, double*, double*, double*, double*, double*, double*, double*, int*> TupleIt;
  typedef thrust::zip_iterator< TupleIt >  ZipIt;
*/
//  int streamID = acc_get_cuda_stream(acc_async_sync);
  /*ZipIt Zend =*/ thrust::remove_if(thrust::cuda::par.on((cudaStream_t)cuda_stream),
                       thrust::make_zip_iterator(thrust::make_tuple(x_star, y_star, z_star,
                                                                    x, y, z, v_x, v_y, v_z)),
                       thrust::make_zip_iterator(thrust::make_tuple(x_star+count, y_star+count, z_star+count,
                                                                    x+count, y+count, z+count, v_x+count, v_y+count, v_z+count)),
                       RemoveOutsideBounds(min, max));

/*  TupleIt endTuple = Zend.get_iterator_tuple();
  printf("num removed: %d\n", x_star + count - thrust::get<0>(endTuple));*/
}
