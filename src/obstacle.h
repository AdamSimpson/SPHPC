#ifndef SPH_SRC_OBSTACLE_H_
#define SPH_SRC_OBSTACLE_H_

// Ehh....use forward declaration everywhere else
// could just store pointers in Obstacle struct
#include "obstacle_cuda.h"
#include "geometry.h"
// Forward Declaration
//struct AABB;
//struct Obstacle_CUDA;

struct Obstacle {
  double origin_x;
  double origin_y;
  double origin_z;

  struct AABB world_bounds; // Bounding box world dimensions

  double SDF_spacing; // Spacing of signed distance field cells

  int SDF_dim_x;
  int SDF_dim_y;
  int SDF_dim_z;

  float *SDF_buffer;

  struct Obstacle_CUDA obstacle_cuda;
};

void AllocInitObstacle(struct Obstacle *obstacle, const char *SDF_file_name);
void FinalizeObstacle(struct Obstacle *obstacle);

#endif
