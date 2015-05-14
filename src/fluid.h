#ifndef SPH_SRC_FLUID_H_
#define SPH_SRC_FLUID_H_

typedef struct FLUID_PARTICLES FluidParticles;

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

struct FLUID_PARTICLES {
  double *x_star;
  double *y_star;
  double *z_star;
  double *x;
  double *y;
  double *z;
  double *v_x;
  double *v_y;
  double *v_z;
  double *dp_x;
  double *dp_y;
  double *dp_z;
  double *w_x;
  double *w_y;
  double *w_z;
  double *density;
  double *lambda;
  int *id; // Id is 'local' index within the fluid particle pointer array
};

////////////////////////////////////////////////
// Function prototypes
////////////////////////////////////////////////
void MoveParticle(FluidParticles *const particles, const int from_index, const int to_index);

void VorticityConfinement(FluidParticles *const fluid_particles,
                          const Params *const params,
                          const Neighbors *const neighbors);

void XSPHViscosity(FluidParticles *const fluid_particles,
                   const Params *const params,
                   const Neighbors *const neighbors);

void ComputeDensities(FluidParticles *const fluid_particles,
                      const Params *const params,
                      const Neighbors *const neighbors);

void ApplyGravity(FluidParticles *const fluid_particles,
                  const Params *const params);

void UpdatePositionStars(FluidParticles *const fluid_particles,
                         const Params *const params,
                         const AABB *const boundary_global);

void UpdatePositions(FluidParticles *const fluid_particles,
                     const Params *const params);

void CalculateLambda(FluidParticles *const fluid_particles,
                     const Params *const params,
                     const Neighbors *const neighbors);

void UpdateDPs(FluidParticles *const fluid_particles,
               const Params *const params,
               const Neighbors *const neighbors);

void IdentifyOOBParticles(FluidParticles *const fluid_particles,
                          Params *const params,
                          Communication *const communication);

void PredictPositions(FluidParticles *const fluid_particles,
                      const Params *const params,
                      const AABB *const boundary_global);

void CheckVelocity(double *const v_x, double *const v_y, double *const v_z);

void UpdateVelocities(FluidParticles *const fluid_particles,
                      const Params *const params);

void BoundaryConditions(FluidParticles *const fluid_particles,
                        const unsigned int i,
                        const AABB *const boundary);

void InitParticles(FluidParticles *const fluid_particles,
                   Params *const params,
                   const AABB *const water);

void AllocateFluid(FluidParticles *particles,
                    const Params *const params);

void FreeFluid(FluidParticles *particles);

#endif
