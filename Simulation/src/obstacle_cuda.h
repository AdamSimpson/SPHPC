#ifndef SPH_SRC_OBSTACLECUDA_H_
#define SPH_SRC_OBSTACLECUDA_H_

#include "cuda_runtime.h"

// Forward struct declarations
struct Obstacle;
struct obstacle_CUDA;
struct particles;

struct Obstacle_CUDA {
  struct cudaArray *SDF_buffer;
  cudaTextureObject_t SDF_tex;
};

// These functions are to be visible to both C and C++
#ifdef __cplusplus
extern "C" {
#endif

void AllocInitObstacle_CUDA(struct Obstacle *obstacle);

void FinalizeObstacle_CUDA(struct Obstacle_CUDA *obstacle_cuda);

void CalculateParticleCollisions(double *restrict x,
                                            double *restrict y,
                                            double *restrict z,
                                            const int particle_count,
                                            const struct Obstacle obstacle);

// These functions are to be visible to both C and C++
#ifdef __cplusplus
}
#endif

#endif
