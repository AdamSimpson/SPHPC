#ifndef SPH_SRC_SIMULATION_H_
#define SPH_SRC_SIMULATION_H_

typedef struct PARAM Params;

#include "geometry.h"

struct PARAM {
  double rest_density;
  double smoothing_radius;
  double g;
  double time_step;
  double k;
  double dq;
  double c;
  double node_start_x; // left x position of node partition
  double node_end_x;   // right x position of node partition
  int number_fluid_particles_global;
  int number_fluid_particles_local;  // Number of non vacant particles
  int max_fluid_particles_local;     // Maximum number for max_fluid_particle_index + halo particles
  int number_halo_particles_left;    // Starting at max_fluid_particle_index
  int number_halo_particles_right;
  int number_steps;
  int rank;
  int num_procs;
}; // Simulation paramaters

void SetParameters(Params *const params,
                   Neighbors *const neighbors,
                   AABB *const boundary_global,
                   AABB *const water_volume_global);

#endif
