#ifndef SPH_SRC_OBSTACLE_H_
#define SPH_SRC_OBSTACLE_H_

#include "obstacle_cuda.h"
#include "geometry.h"

struct Obstacle {
  char SDF_file_name[128];
  double origin_x;
  double origin_y;
  double origin_z;

  double model_length;
  double model_width;
  double model_height;

  struct AABB world_bounds; // Bounding box world dimensions

  double SDF_spacing; // Spacing of signed distance field cells

  int SDF_dim_x;
  int SDF_dim_y;
  int SDF_dim_z;

  float *SDF_buffer;

  struct Obstacle_CUDA obstacle_cuda;
};

void AllocInitObstacle(struct Obstacle *obstacle);
void FinalizeObstacle(struct Obstacle *obstacle);

#endif
