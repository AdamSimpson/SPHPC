#ifndef SPH_SRC_THRUST_C_H_
#define SPH_SRC_THRUST_C_H_

// These functions are to be visible to both C and C++
#ifdef __cplusplus
extern "C" {
#endif

void CopyIfLessThan(const double min,
                    const int *const input,
                    const int input_count,
                    const double *const stencil,
                    int *const output,
                    int *const output_count,
                    const void *const cuda_stream);

void CopyIfLessThanOrEqual(const double min,
                           const int *const input,
                           const int input_count,
                           const double *const stencil,
                           int *const output,
                           int *const output_count,
                           const void *const cuda_stream);

void CopyIfGreaterThan(const double max,
                       const int *const input,
                       const int input_count,
                       const double *const stencil,
                       int *const output,
                       int *const output_count,
                       const void *const cuda_stream);

void CopyIfGreaterThanOrEqual(const double max,
                       const int *const input,
                       const int input_count,
                       const double *const stencil,
                       int *const output,
                       int *const output_count,
                       const void *const cuda_stream);

void RemoveIfOutsideBounds(const double min, const double max,
                           int *const input,
                           const int input_count,
                           const double *const stencil,
                           const void *const cuda_stream);

void SortByKey(unsigned int *const keys,
               unsigned int *const values,
               const int count,
               const void *const cuda_stream);

void FindLowerBounds(const unsigned int *const hash_values,
                                const int value_count,
                                const int search_count,
                                unsigned int *const end_indices,
                                const void *const cuda_stream);

void SortParticlesByKey(unsigned int *const keys,
                                 const int count,
                                 double *x_star, double *y_star, double *z_star,
                                 double *x, double *y, double *z,
                                 double *v_x, double *v_y, double *v_z,
                                 const void *const cuda_stream);

void FindUpperBounds(const unsigned int *const hash_values,
                                const int value_count,
                                const int search_count,
                                unsigned int *const start_indices,
                                const void *const cuda_stream);

void RemoveIfOutsideBounds2(const double min, const double max,
                                      int *const input,
                                      const int input_count,
                                      const void *const cuda_stream);

void RemoveIfOutsideBounds3(const double min, const double max,
                                      const int count,
                                      double *x_star, double *y_star, double *z_star,
                                      double *x, double *y, double *z,
                                      double *v_x, double *v_y, double *v_z,
                                      const void *const cuda_stream);

#ifdef __cplusplus
}
#endif

#endif
