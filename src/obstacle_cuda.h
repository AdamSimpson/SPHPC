#ifndef SPH_SRC_OBSTACLECUDA_H_
#define SPH_SRC_OBSTACLECUDA_H_

#include "cuda_runtime.h"

// Forward struct declarations
struct Obstacle;

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

bool InsideObstacle(const struct Obstacle *obstacle,
                    const double x, const double y, const double z);

// These functions are to be visible to both C and C++
#ifdef __cplusplus
}
#endif

#endif
