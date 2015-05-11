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
void VorticityConfinement(FluidParticle *const fluid_particles,
                          const Params *const params,
                          const Neighbors *const neighbors);

void XSPHViscosity(FluidParticle *const fluid_particles,
                   const Params *const params,
                   const Neighbors *const neighbors);

void ComputeDensities(FluidParticle *const fluid_particles,
                      const Params *const params,
                      const Neighbors *const neighbors);

void ApplyGravity(FluidParticle *const fluid_particles,
                  const Params *const params);

void UpdatePositionStars(FluidParticle *const fluid_particles,
                         const Params *const params,
                         const AABB *const boundary_global);

void UpdatePositions(FluidParticle *const fluid_particles,
                     const Params *const params);

void CalculateLambda(FluidParticle *const fluid_particles,
                     const Params *const params,
                     const Neighbors *const neighbors);

void UpdateDPs(FluidParticle *const fluid_particles,
               const Params *const params,
               const Neighbors *const neighbors);

void IdentifyOOBParticles(FluidParticle *const fluid_particles,
                          Params *const params,
                          Communication *const communication);

void PredictPositions(FluidParticle *const fluid_particles,
                      const Params *const params,
                      const AABB *const boundary_global);

void CheckVelocity(double *const v_x, double *const v_y, double *const v_z);

void UpdateVelocities(FluidParticle *const fluid_particles,
                      const Params *const params);

void BoundaryConditions(FluidParticle *const fluid_particles,
                        const unsigned int i,
                        const AABB *const boundary);

void InitParticles(FluidParticle *const fluid_particles,
                   Params *const params,
                   const AABB *const water);

void AllocateFluid(FluidParticle **fluid_particles,
                   const Params *const params);

#endif
