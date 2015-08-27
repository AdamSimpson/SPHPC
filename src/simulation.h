#ifndef SPH_SRC_SIMULATION_H_
#define SPH_SRC_SIMULATION_H_

// Forward Declarations
struct Particles;
struct Neighbors;
struct AABB;
struct Obstacle;

struct Params {
  double smoothing_radius;
  double W_norm;
  double DelW_norm;
  double g;
  double time_step;
  double k;
  double dq;
  double c;
  double node_start_x; // left x position of node partition
  double node_end_x;   // right x position of node partition
  int number_steps;
  int rank;
  int proc_count;
}; // Simulation paramaters

void SetParameters(struct Params *const params,
                   struct Particles *const particles,
                   struct Neighbors *const neighbors,
                   struct AABB *const boundary_global,
                   struct AABB *const water_volume_global,
                   struct Obstacle *const obstacle);

#endif
