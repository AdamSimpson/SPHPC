#ifndef SPH_SRC_FLUID_H_
#define SPH_SRC_FLUID_H_

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Forward Declaration
struct Params;
struct Neighbors;
struct AABB;
struct Communication;

////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////

struct FluidParticles {
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
void MoveParticle(struct FluidParticles *const particles, const int from_index, const int to_index);

void VorticityConfinement(struct FluidParticles *const fluid_particles,
                          const struct Params *const params,
                          const struct Neighbors *const neighbors);

void XSPHViscosity(struct FluidParticles *const fluid_particles,
                   const struct Params *const params,
                   const struct Neighbors *const neighbors);

void ComputeDensities(struct FluidParticles *const fluid_particles,
                      const struct Params *const params,
                      const struct Neighbors *const neighbors);

void ApplyGravity(struct FluidParticles *const fluid_particles,
                  const struct Params *const params);

void UpdatePositionStars(struct FluidParticles *const fluid_particles,
                         const struct Params *const params,
                         const struct AABB *const boundary_global);

void UpdatePositions(struct FluidParticles *const fluid_particles,
                     const struct Params *const params);

void CalculateLambda(struct FluidParticles *const fluid_particles,
                     const struct Params *const params,
                     const struct Neighbors *const neighbors);

void UpdateDPs(struct FluidParticles *const fluid_particles,
               const struct Params *const params,
               const struct Neighbors *const neighbors);

void IdentifyOOBParticles(struct FluidParticles *const fluid_particles,
                          struct Params *const params,
                          struct Communication *const communication);

void PredictPositions(struct FluidParticles *const fluid_particles,
                      const struct Params *const params,
                      const struct AABB *const boundary_global);

void CheckVelocity(double *const v_x, double *const v_y, double *const v_z);

void UpdateVelocities(struct FluidParticles *const fluid_particles,
                      const struct Params *const params);

void BoundaryConditions(struct FluidParticles *const fluid_particles,
                        const unsigned int i,
                        const struct AABB *const boundary);

void InitParticles(struct FluidParticles *const fluid_particles,
                   struct Params *const params,
                   const struct AABB *const water);

void AllocateFluid(struct FluidParticles *particles,
                  const struct Params *const params);

void FreeFluid(struct FluidParticles *particles);

#endif
