#ifndef SPH_SRC_OUT_OF_BOUNDS_H_
#define SPH_SRC_OUT_OF_BOUNDS_H_

// Forward Declarations
struct Particles;
struct Communication;
struct Params;

// These functions are to be visible to both C and C++
#ifdef __cplusplus
extern "C" {
#endif

void PackOOB(struct Params *const params,
             struct Particles *const particles,
             struct Communication *const communication);

void CopyIfLessThan(const double min,
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

void RemoveIfOutsideBounds(const double min, const double max,
                           int *const input,
                           const int input_count,
                           const double *const stencil);

#ifdef __cplusplus
}
#endif

#endif
