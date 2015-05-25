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

#ifdef __cplusplus
}
#endif

#endif
