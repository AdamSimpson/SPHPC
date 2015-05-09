#ifndef SPH_SRC_FLUID_H_
#define SPH_SRC_FLUID_H_

typedef struct FLUID_PARTICLE fluid_particle_t;

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "neighbors.h"
#include "fileio.h"
#include "geometry.h"
#include "communication.h"
#include "simulation.h"

////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////

struct FLUID_PARTICLE {
    double x_star;
    double y_star;
    double z_star;
    double x;
    double y;
    double z;
    double v_x;
    double v_y;
    double v_z;
    double dp_x;
    double dp_y;
    double dp_z;
    double w_x;
    double w_y;
    double w_z;
    double density;
    double lambda;
    int id; // Id is 'local' index within the fluid particle pointer array
};

////////////////////////////////////////////////
// Function prototypes
////////////////////////////////////////////////
void vorticity_confinement(fluid_particle_t *const fluid_particles,
                           const param_t *const params,
                           const neighbors_t *const neighbors);

void XSPH_viscosity(fluid_particle_t *const fluid_particles,
                    const param_t *const params,
                    const neighbors_t *const neighbors);

void compute_densities(fluid_particle_t *const fluid_particles,
                       const param_t *const params,
                       const neighbors_t *const neighbors);

void apply_gravity(fluid_particle_t *const fluid_particles,
                  const param_t *const params);

void update_dp_positions(fluid_particle_t *const fluid_particles,
                         const param_t *const params,
                         const AABB_t *const boundary_global);

void update_positions(fluid_particle_t *const fluid_particles,
                      const param_t *const params);

void calculate_lambda(fluid_particle_t *const fluid_particles,
                      const param_t *const params,
                      const neighbors_t *const neighbors);

void update_dp(fluid_particle_t *const fluid_particles,
               const param_t *const params,
               const neighbors_t *const neighbors);

void identify_oob_particles(fluid_particle_t *const fluid_particles,
                            param_t *const params,
                            communication_t *const communication);

void predict_positions(fluid_particle_t *const fluid_particles,
                       const param_t *const params,
                       const AABB_t *const boundary_global);

void check_velocity(double *const v_x, double *const v_y, double *const v_z);

void update_velocities(fluid_particle_t *const fluid_particles,
                       const param_t *const params);

void boundary_conditions(fluid_particle_t *const fluid_particles,
                         const unsigned int i,
                         const AABB_t *const boundary);

void init_particles(fluid_particle_t *const fluid_particles,
                    param_t *const params,
                    const AABB_t *const water);

void allocate_fluid(fluid_particle_t **fluid_particles,
                    const param_t *const params);

#endif
