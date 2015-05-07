#ifndef SPH_SRC_SIMULATION_H_
#define SPH_SRC_SIMULATION_H_

typedef struct PARAM param_t;

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
    int nprocs;
}; // Simulation paramaters

void set_parameters(param_t *params, neighbors_t *neighbors,
                    AABB_t *boundary_global, AABB_t* water_volume_global);

#endif
