#ifndef SPH_SRC_FLUID_H_
#define SPH_SRC_FLUID_H_

// Forward Declaration
struct Params;
struct Neighbors;
struct AABB;
struct Communication;

////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////

struct Particles {
  double *restrict x_star;
  double *restrict y_star;
  double *restrict z_star;
  double *restrict x;
  double *restrict y;
  double *restrict z;
  double *restrict v_x;
  double *restrict v_y;
  double *restrict v_z;
  double *restrict dp_x;
  double *restrict dp_y;
  double *restrict dp_z;
  double *restrict w_x;
  double *restrict w_y;
  double *restrict w_z;
  double *restrict density;
  double *restrict lambda;
  int    *restrict id; // Id is 'local' index within the fluid particle pointer array

  int global_count;
  int local_count;
  int max_local;  // Maximum number for max_fluid_particle_index + halo particles
  int halo_count_left;
  int halo_count_right;
};

////////////////////////////////////////////////
// Function prototypes
////////////////////////////////////////////////
void MoveParticle(struct Particles *const particles, const int from_index, const int to_index);

void VorticityConfinement(struct Particles *const fluid_particles,
                          const struct Params *const params,
                          const struct Neighbors *const neighbors);

void XSPHViscosity(struct Particles *const fluid_particles,
                   const struct Params *const params,
                   const struct Neighbors *const neighbors);

void ComputeDensities(struct Particles *const fluid_particles,
                      const struct Params *const params,
                      const struct Neighbors *const neighbors);

void ApplyGravity(struct Particles *const fluid_particles,
                  const struct Params *const params);

void UpdatePositionStars(struct Particles *const fluid_particles,
                         const struct AABB *const boundary_global);

void UpdatePositions(struct Particles *const fluid_particles);

void CalculateLambda(struct Particles *const fluid_particles,
                     const struct Params *const params,
                     const struct Neighbors *const neighbors);

void UpdateDPs(struct Particles *const fluid_particles,
               const struct Params *const params,
               const struct Neighbors *const neighbors);

void IdentifyOOBParticles(struct Particles *const fluid_particles,
                          struct Params *const params,
                          struct Communication *const communication);

void PredictPositions(struct Particles *const fluid_particles,
                      const struct Params *const params,
                      const struct AABB *const boundary_global);

void CheckVelocity(double *const v_x, double *const v_y, double *const v_z);

void UpdateVelocities(struct Particles *const fluid_particles,
                      const struct Params *const params);

void BoundaryConditions(struct Particles *const fluid_particles,
                        const unsigned int i,
                        const struct AABB *const boundary);

void InitParticles(struct Particles *const fluid_particles,
                   struct Params *const params,
                   const struct AABB *const water);

void AllocateFluid(struct Particles *particles);

void FreeFluid(struct Particles *particles);

#endif
