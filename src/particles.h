#ifndef SPH_SRC_PARTICLES_H_
#define SPH_SRC_PARTICLES_H_

// Forward Declaration
struct Params;
struct Neighbors;
struct AABB;
struct Communication;
struct Obstacle;

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

  double rest_density;
  double rest_mass;
  int global_count; // Global number of particles in simulation
  int local_count; // Particles within node bounds, excludes halo particles
  int max_local;  // Maximum number of local and halo particles
  int halo_count_left;
  int halo_count_right;
};

void ComputeVorticity(struct Particles *restrict particles,
                      const struct Params *restrict params,
                      const struct Neighbors *restrict neighbors);

void ApplyVorticityConfinement(struct Particles *restrict fluid_particles,
                               const struct Params *restrict params,
                               const struct Neighbors *restrict neighbors);

void ApplyViscosity(struct Particles *restrict fluid_particles,
                    const struct Params *restrict params,
                    const struct Neighbors *restrict neighbors);

void ComputeDensities(struct Particles *restrict fluid_particles,
                      const struct Params *restrict params,
                      const struct Neighbors *restrict neighbors);

void ApplyGravity(struct Particles *restrict fluid_particles,
                  const struct Params *restrict params);

void UpdatePositionStars(struct Particles *restrict particles,
                         const struct AABB *restrict boundary_global,
                         const struct Obstacle *restrict obstacle,
                         const struct Params *params);

void UpdatePositions(struct Particles *restrict fluid_particles);

void ComputeLambda(struct Particles *restrict fluid_particles,
                     const struct Params *restrict params,
                     const struct Neighbors *restrict neighbors);

void UpdateDPs(struct Particles *restrict fluid_particles,
               const struct Params *restrict params,
               const struct Neighbors *restrict neighbors);

void PredictPositions(struct Particles *restrict particles,
                      const struct Params *restrict params,
                      const struct AABB *restrict boundary_global);

void UpdateVelocities(struct Particles *restrict fluid_particles,
                      const struct Params *restrict params);

void ApplyBoundaryConditions(double *restrict x, double *restrict y, double *restrictz ,
                             const struct AABB *restrict boundary,
                             const struct Params *params);

void PrintAverageDensity(struct Particles *restrict particles);

void AllocInitParticles(struct Particles *restrict particles,
                        struct Params *restrict params,
                        struct AABB *restrict fluid_volume_initial);

void FinalizeParticles(struct Particles *restrict particles);

#endif
