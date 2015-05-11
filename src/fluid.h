#ifndef SPH_SRC_FLUID_H_
#define SPH_SRC_FLUID_H_

typedef struct FLUID_PARTICLE FluidParticle;

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
void vorticity_confinement(FluidParticle *const fluid_particles,
                           const Params *const params,
                           const Neighbors *const neighbors);

void XSPH_viscosity(FluidParticle *const fluid_particles,
                    const Params *const params,
                    const Neighbors *const neighbors);

void compute_densities(FluidParticle *const fluid_particles,
                       const Params *const params,
                       const Neighbors *const neighbors);

void apply_gravity(FluidParticle *const fluid_particles,
                  const Params *const params);

void update_dp_positions(FluidParticle *const fluid_particles,
                         const Params *const params,
                         const AABB *const boundary_global);

void update_positions(FluidParticle *const fluid_particles,
                      const Params *const params);

void calculate_lambda(FluidParticle *const fluid_particles,
                      const Params *const params,
                      const Neighbors *const neighbors);

void update_dp(FluidParticle *const fluid_particles,
               const Params *const params,
               const Neighbors *const neighbors);

void identify_oob_particles(FluidParticle *const fluid_particles,
                            Params *const params,
                            Communication *const communication);

void predict_positions(FluidParticle *const fluid_particles,
                       const Params *const params,
                       const AABB *const boundary_global);

void check_velocity(double *const v_x, double *const v_y, double *const v_z);

void update_velocities(FluidParticle *const fluid_particles,
                       const Params *const params);

void boundary_conditions(FluidParticle *const fluid_particles,
                         const unsigned int i,
                         const AABB *const boundary);

void init_particles(FluidParticle *const fluid_particles,
                    Params *const params,
                    const AABB *const water);

void allocate_fluid(FluidParticle **fluid_particles,
                    const Params *const params);

#endif
