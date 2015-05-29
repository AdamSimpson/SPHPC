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
                    int *const output_count);

void CopyIfLessThanOrEqual(const double min,
                           const int *const input,
                           const int input_count,
                           const double *const stencil,
                           int *const output,
                           int *const output_count);

void CopyIfGreaterThan(const double max,
                       const int *const input,
                       const int input_count,
                       const double *const stencil,
                       int *const output,
                       int *const output_count);

void CopyIfGreaterThanOrEqual(const double max,
                       const int *const input,
                       const int input_count,
                       const double *const stencil,
                       int *const output,
                       int *const output_count);

void RemoveIfOutsideBounds(const double min, const double max,
                           int *const input,
                           const int input_count,
                           const double *const stencil);

void SortByKey(unsigned int *const keys,
               unsigned int *const values,
               const int count);

#ifdef __cplusplus
}
#endif

#endif
