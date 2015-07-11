#ifndef SPH_SRC_PARTICLES_H_
#define SPH_SRC_PARTICLES_H_

// Forward Declaration
struct Params;
struct Neighbors;
struct AABB;
struct Communication;

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
  int    *restrict id;

  int global_count; // Global number of particles in simulation
  int local_count; // Particles within node bounds, excludes halo particles
  int max_local;  // Maximum number of local and halo particles
  int halo_count_left;
  int halo_count_right;
};

void ComputeVorticity(struct Particles *const particles,
                               const struct Params *const params,
                               const struct Neighbors *const neighbors);

void ApplyVorticityConfinement(struct Particles *const fluid_particles,
                          const struct Params *const params,
                          const struct Neighbors *const neighbors);

void ApplyViscosity(struct Particles *const fluid_particles,
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

void ComputeLambda(struct Particles *const fluid_particles,
                     const struct Params *const params,
                     const struct Neighbors *const neighbors);

void UpdateDPs(struct Particles *const fluid_particles,
               const struct Params *const params,
               const struct Neighbors *const neighbors);

void PredictPositions(struct Particles *const fluid_particles,
                      const struct Params *const params,
                      const struct AABB *const boundary_global);

void UpdateVelocities(struct Particles *const fluid_particles,
                      const struct Params *const params);

void ApplyBoundaryConditions(double *x, double *y, double*z,
                        const struct AABB *const boundary);

void PrintAverageDensity(struct Particles *particles);

void AllocInitParticles(struct Particles *particles,
                        struct Params *params,
                        struct AABB *fluid_volume_initial);

void FinalizeParticles(struct Particles *particles);

#endif
