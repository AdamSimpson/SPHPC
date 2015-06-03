#include "particles.h"
#include "simulation.h"
#include "communication.h"
#include "neighbors.h"
#include "geometry.h"
#include "debug.h"
#include "safe_alloc.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

///////////////////////////////////////////////////////////////////////////
// Smoothing Kernels
// Don't use pow as it's POWerfully slow
///////////////////////////////////////////////////////////////////////////

// (h^2 - r^2)^3 normalized in 3D (poly6)
static inline double W(const double r, const double h) {
  if(r > h)
    return 0.0;

  const double C = 315.0/(64.0*M_PI* h*h*h*h*h*h*h*h*h);
  const double W = C*(h*h-r*r)*(h*h-r*r)*(h*h-r*r);
  return W;
}

// Gradient (h-r)^3 normalized in 3D (Spikey) magnitude
// Need to multiply by r/|r|
static inline double DelW(const double r, const double h) {
  if(r > h)
    return 0.0;

  const double C = -45.0/(M_PI * h*h*h*h*h*h);
  const double DelW = C*(h-r)*(h-r);
  return DelW;
}

void CopyParticle(struct Particles *const particles,
                  const int from_index, const int to_index) {
  particles->x_star[to_index]  = particles->x_star[from_index];
  particles->y_star[to_index]  = particles->y_star[from_index];
  particles->z_star[to_index]  = particles->z_star[from_index];
  particles->x[to_index]       = particles->x[from_index];
  particles->y[to_index]       = particles->y[from_index];
  particles->z[to_index]       = particles->z[from_index];
  particles->v_x[to_index]     = particles->v_x[from_index];
  particles->v_y[to_index]     = particles->v_y[from_index];
  particles->v_z[to_index]     = particles->v_z[from_index];
  particles->dp_x[to_index]    = particles->dp_x[from_index];
  particles->dp_y[to_index]    = particles->dp_y[from_index];
  particles->dp_z[to_index]    = particles->dp_z[from_index];
  particles->w_x[to_index]     = particles->w_x[from_index];
  particles->w_y[to_index]     = particles->w_y[from_index];
  particles->w_z[to_index]     = particles->w_z[from_index];
  particles->density[to_index] = particles->density[from_index];
  particles->lambda[to_index]  = particles->lambda[from_index];
  particles->id[to_index]      = to_index;
}

////////////////////////////////////////////////////////////////////////////
// Particle attribute computations
////////////////////////////////////////////////////////////////////////////

void ApplyVorticityConfinement(struct Particles *const particles,
                           const struct Params *const params,
                           const struct Neighbors *const neighbors) {

  const double dt = params->time_step;
  const double eps = 5.0;

  // Calculate vorticy at each particle
  for (int i=0; i<particles->local_count; ++i) {
    const int p_index = i;
    const struct NeighborBucket *const n = &neighbors->neighbor_buckets[i];

    double vort_x = 0.0;
    double vort_y = 0.0;
    double vort_z = 0.0;

    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];
      const double x_diff = particles->x_star[p_index]
                          - particles->x_star[q_index];
      const double y_diff = particles->y_star[p_index]
                          - particles->y_star[q_index];
      const double z_diff = particles->z_star[p_index]
                          - particles->z_star[q_index];

      const double vx_diff = particles->v_x[q_index] - particles->v_x[p_index];
      const double vy_diff = particles->v_y[q_index] - particles->v_y[p_index];
      const double vz_diff = particles->v_z[q_index] - particles->v_z[p_index];

      double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
      if(r_mag < 0.0001)
        r_mag = 0.0001;

        const double dw = DelW(r_mag, params->smoothing_radius);

        const double dw_x = dw*x_diff/r_mag;
        const double dw_y = dw*y_diff/r_mag;
        const double dw_z = dw*z_diff/r_mag;

        vort_x +=  vy_diff*dw_z - vz_diff*dw_y;
        vort_y +=  vz_diff*dw_x - vx_diff*dw_z;
        vort_z +=  vx_diff*dw_y - vy_diff*dw_x;
    }

    particles->w_x[p_index] = vort_x;
    particles->w_y[p_index] = vort_y;
    particles->w_z[p_index] = vort_z;
  }

  // Apply vorticity confinement
  for (int i=0; i<particles->local_count; ++i) {
    const int p_index = i;
    const struct NeighborBucket *const n = &neighbors->neighbor_buckets[i];

    double eta_x  = 0.0;
    double eta_y  = 0.0;
    double eta_z  = 0.0;

    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];

      const double x_diff = particles->x_star[p_index]
                          - particles->x_star[q_index];
      const double y_diff = particles->y_star[p_index]
                          - particles->y_star[q_index];
      const double z_diff = particles->z_star[p_index]
                          - particles->z_star[q_index];

      double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
      if(r_mag < 0.0001)
        r_mag = 0.0001;

      const double dw = DelW(r_mag, params->smoothing_radius);

      const double dw_x = dw*x_diff/r_mag;
      const double dw_y = dw*y_diff/r_mag;
      const double dw_z = dw*z_diff/r_mag;

      const double vort_mag = sqrt(particles->w_x[q_index]*particles->w_x[q_index]
                                 + particles->w_y[q_index]*particles->w_y[q_index]
                                 + particles->w_z[q_index]*particles->w_z[q_index]);

      eta_x += vort_mag*dw_x;
      eta_y += vort_mag*dw_y;
      eta_z += vort_mag*dw_z;
    }

    double eta_mag = sqrt(eta_x*eta_x + eta_y*eta_y + eta_z*eta_z);
    if(eta_mag < 0.0001)
      eta_mag = 0.0001;

    const double N_x = eta_x/eta_mag;
    const double N_y = eta_y/eta_mag;
    const double N_z = eta_z/eta_mag;

    const double w_x = particles->w_x[p_index];
    const double w_y = particles->w_y[p_index];
    const double w_z = particles->w_z[p_index];

    particles->v_x[p_index] += dt*eps * (N_y*w_z - N_z*w_y);
    particles->v_y[p_index] += dt*eps * (N_z*w_x - N_x*w_z);
    particles->v_z[p_index] += dt*eps * (N_x*w_y - N_y*w_x);
  }
}

void ApplyViscosity(struct Particles *const particles,
                   const struct Params *const params,
                   const struct Neighbors *const neighbors) {

  const double c = params->c;
  const double h = params->smoothing_radius;

  for (int i=0; i < particles->local_count; ++i) {
    const int p_index = i;
    const struct NeighborBucket *const n = &neighbors->neighbor_buckets[i];

    double partial_sum_x = 0.0;
    double partial_sum_y = 0.0;
    double partial_sum_z = 0.0;

    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];
      const double x_diff = particles->x_star[p_index]
                          - particles->x_star[q_index];
      const double y_diff = particles->y_star[p_index]
                          - particles->y_star[q_index];
      const double z_diff = particles->z_star[p_index]
                          - particles->z_star[q_index];

      const double vx_diff = particles->v_x[q_index]
                           - particles->v_x[p_index];
      const double vy_diff = particles->v_y[q_index]
                           - particles->v_y[p_index];
      const double vz_diff = particles->v_z[q_index]
                           - particles->v_z[p_index];

      const double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff * z_diff*z_diff);

      const double w = W(r_mag, h);

      // http://mmacklin.com/pbf_sig_preprint.pdf is missing 1/sigma contribution
      // see: http://www.cs.ubc.ca/~rbridson/docs/schechter-siggraph2012-ghostsph.pdf
      partial_sum_x += vx_diff * w / particles->density[q_index];
      partial_sum_y += vy_diff * w / particles->density[q_index];
      partial_sum_z += vz_diff * w / particles->density[q_index];
    }

    partial_sum_x *= c;
    partial_sum_y *= c;
    partial_sum_z *= c;

    particles->v_x[p_index] += partial_sum_x;
    particles->v_y[p_index] += partial_sum_y;
    particles->v_z[p_index] += partial_sum_z;
  }
}

void ComputeDensities(struct Particles *const particles,
                      const struct Params *const params,
                      const struct Neighbors *const neighbors) {

  const double h = params->smoothing_radius;

  for (int i=0; i<particles->local_count; ++i) {
    const int p_index = i;
    const struct NeighborBucket *const n = &neighbors->neighbor_buckets[i];

    // Own contribution to density
    double density = W(0.0, h);

    // Neighbor contribution
    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];
      const double x_diff = particles->x_star[p_index]
                          - particles->x_star[q_index];
      const double y_diff = particles->y_star[p_index]
                          - particles->y_star[q_index];
      const double z_diff = particles->z_star[p_index]
                          - particles->z_star[q_index];

      const double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
      density += W(r_mag, h);
    }

    // Update particle density
    particles->density[p_index] = density;
  }
}

void ApplyGravity(struct Particles *const particles,
                  const struct Params *const params) {
  const double dt = params->time_step;
  const double g = -params->g;

  for (int i=0; i<(particles->local_count); ++i) {
    particles->v_y[i] += g*dt;
  }
}

void UpdatePositionStars(struct Particles *const particles,
                         const struct AABB *const boundary_global) {

  for (int i=0; i<(particles->local_count); ++i) {
    particles->x_star[i] += particles->dp_x[i];
    particles->y_star[i] += particles->dp_y[i];
    particles->z_star[i] += particles->dp_z[i];

    // Enforce boundary conditions
    ApplyBoundaryConditions(particles, i, boundary_global);
  }

}

void UpdatePositions(struct Particles *const particles) {

  for (int i=0; i<particles->local_count; ++i) {
    particles->x[i] = particles->x_star[i];
    particles->y[i] = particles->y_star[i];
    particles->z[i] = particles->z_star[i];
  }

}

void ComputeLambda(struct Particles *const particles,
                     const struct Params *const params,
                     const struct Neighbors *const neighbors) {

  for (int i=0; i<particles->local_count; ++i) {
    const int p_index = i;
    const struct NeighborBucket *const n = &neighbors->neighbor_buckets[i];
    const double Ci = particles->density[p_index]
                     /params->rest_density - 1.0;
    double sum_C = 0.0;
    double sum_grad_x = 0.0;
    double sum_grad_y = 0.0;
    double sum_grad_z = 0.0;

    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];
      const double x_diff = particles->x_star[p_index]
                          - particles->x_star[q_index];
      const double y_diff = particles->y_star[p_index]
                          - particles->y_star[q_index];
      const double z_diff = particles->z_star[p_index]
                          - particles->z_star[q_index];

      double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
      if (r_mag < 0.0001)
        r_mag = 0.0001;

      const double grad = DelW(r_mag, params->smoothing_radius);

      const double grad_x = grad*x_diff/r_mag;
      const double grad_y = grad*y_diff/r_mag;
      const double grad_z = grad*z_diff/r_mag;
      sum_grad_x += grad_x;
      sum_grad_y += grad_y;
      sum_grad_z += grad_z;
      // Add k = j contribution
      sum_C += (grad_x*grad_x + grad_y*grad_y + grad_z*grad_z);
    }

    // k = i contribution
    sum_C += sum_grad_x*sum_grad_x
           + sum_grad_y*sum_grad_y
           + sum_grad_z*sum_grad_z;

    sum_C *= (1.0/(params->rest_density*params->rest_density));

    const double epsilon = 1.0;
    particles->lambda[i] = -Ci/(sum_C + epsilon);
  }
}

void UpdateDPs(struct Particles *const particles,
               const struct Params *const params,
               const struct Neighbors *const neighbors) {

  const double h = params->smoothing_radius;
  const double k = params->k;
  const double dq = params->dq;
  const double Wdq = W(dq, h);

  for (int i=0; i<particles->local_count; ++i) {
    const int p_index = i;
    const struct NeighborBucket *const n = &neighbors->neighbor_buckets[i];

    double dp_x = 0.0;
    double dp_y = 0.0;
    double dp_z = 0.0;

    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];
      const double x_diff = particles->x_star[p_index]
                          - particles->x_star[q_index];
      const double y_diff = particles->y_star[p_index]
                          - particles->y_star[q_index];
      const double z_diff = particles->z_star[p_index]
                          - particles->z_star[q_index];

      double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
      if(r_mag < 0.0001)
        r_mag = 0.0001;

      const double WdWdq = W(r_mag, h)/Wdq;
      const double s_corr = -k*WdWdq*WdWdq*WdWdq*WdWdq;
      const double dp = (particles->lambda[p_index]
                       + particles->lambda[q_index] + s_corr)
                       * DelW(r_mag, h);

      dp_x += dp * x_diff/r_mag;
      dp_y += dp * y_diff/r_mag;
      dp_z += dp * z_diff/r_mag;
    }

    particles->dp_x[p_index] = dp_x/params->rest_density;
    particles->dp_y[p_index] = dp_y/params->rest_density;
    particles->dp_z[p_index] = dp_z/params->rest_density;
  }

}

void PredictPositions(struct Particles *const particles,
                      const struct Params *const params,
                      const struct AABB *const boundary_global) {
  const double dt = params->time_step;

  for (int i=0; i<particles->local_count; ++i) {
    particles->x_star[i] = particles->x[i]
                         +(particles->v_x[i] * dt);
    particles->y_star[i] = particles->y[i]
                         +(particles->v_y[i] * dt);
    particles->z_star[i] = particles->z[i]
                         +(particles->v_z[i] * dt);

    // Enforce boundary conditions before hash
    // Otherwise predicted position can blow up hash
    ApplyBoundaryConditions(particles, i, boundary_global);
  }
}

void CheckVelocity(double *const v_x, double *const v_y, double *const v_z) {
  const double v_max = 30.0;

  if (*v_x > v_max)
    *v_x = v_max;
  else if (*v_x < -v_max)
    *v_x = -v_max;

  if (*v_y > v_max)
    *v_y = v_max;
  else if (*v_y < -v_max)
    *v_y = -v_max;

  if (*v_z > v_max)
    *v_z = v_max;
  else if (*v_z < -v_max)
    *v_z = -v_max;

}

// Update particle position and check boundary
void UpdateVelocities(struct Particles *const particles,
                      const struct Params *const params) {

  const double dt = params->time_step;

  // Update local and halo particles, update halo so that XSPH visc. is correct
  int total_particles = particles->local_count
                      + particles->halo_count_left
                      + particles->halo_count_right;
  for(int i=0; i<total_particles; ++i) {
    double v_x = (particles->x_star[i] - particles->x[i])/dt;
    double v_y = (particles->y_star[i] - particles->y[i])/dt;
    double v_z = (particles->z_star[i] - particles->z[i])/dt;

    CheckVelocity(&v_x, &v_y, &v_z);

    particles->v_x[i] = v_x;
    particles->v_y[i] = v_y;
    particles->v_z[i] = v_z;
  }
}

// Assume AABB with min point being axis origin
void ApplyBoundaryConditions(struct Particles *const particles,
                        const unsigned int i,
                        const struct AABB *const boundary) {

  // Make sure object is not outside boundary
  // The particle must not be equal to boundary max or
  // hash potentially won't pick it up
  // as the particle will in the 'next' after last bin
  if (particles->x_star[i]  < boundary->min_x)
    particles->x_star[i] = boundary->min_x;
  else if(particles->x_star[i]  > boundary->max_x)
    particles->x_star[i] = boundary->max_x-0.00001;

  if(particles->y_star[i]  <  boundary->min_y)
    particles->y_star[i] = boundary->min_y;
  else if(particles->y_star[i]  > boundary->max_y)
    particles->y_star[i] = boundary->max_y-0.00001;

  if(particles->z_star[i]  <  boundary->min_z)
    particles->z_star[i] = boundary->min_z;
  else if(particles->z_star[i]  > boundary->max_z)
    particles->z_star[i] = boundary->max_z-0.00001;

}

void AllocInitParticles(struct Particles *particles,
                        struct Params *params,
                        struct AABB *fluid_volume_initial) {

  const size_t num_particles = particles->max_local;

  particles->x_star = SAFE_ALLOC(num_particles, sizeof(double));
  particles->y_star = SAFE_ALLOC(num_particles, sizeof(double));
  particles->z_star = SAFE_ALLOC(num_particles, sizeof(double));

  particles->x = SAFE_ALLOC(num_particles, sizeof(double));
  particles->y = SAFE_ALLOC(num_particles, sizeof(double));
  particles->z = SAFE_ALLOC(num_particles, sizeof(double));

  particles->v_x = SAFE_ALLOC(num_particles, sizeof(double));
  particles->v_y = SAFE_ALLOC(num_particles, sizeof(double));
  particles->v_z = SAFE_ALLOC(num_particles, sizeof(double));

  particles->dp_x = SAFE_ALLOC(num_particles, sizeof(double));
  particles->dp_y = SAFE_ALLOC(num_particles, sizeof(double));
  particles->dp_z = SAFE_ALLOC(num_particles, sizeof(double));

  particles->w_x = SAFE_ALLOC(num_particles, sizeof(double));
  particles->w_y = SAFE_ALLOC(num_particles, sizeof(double));
  particles->w_z = SAFE_ALLOC(num_particles, sizeof(double));


  particles->density = SAFE_ALLOC(num_particles, sizeof(double));
  particles->lambda = SAFE_ALLOC(num_particles, sizeof(double));

  particles->id = SAFE_ALLOC(num_particles, sizeof(int));

  ConstructFluidVolume(particles, params, fluid_volume_initial);

  // Initialize particle values
  for (int i=0; i<particles->local_count; ++i) {
    particles->x_star[i] = particles->x[i];
    particles->y_star[i] = particles->y[i];
    particles->z_star[i] = particles->z[i];
    particles->density[i] = params->rest_density;
  }

  #pragma acc enter data copyin(           \
    particles->x_star[0:num_particles],  \
    particles->y_star[0:num_particles],  \
    particles->z_star[0:num_particles],  \
    particles->x[0:num_particles],       \
    particles->y[0:num_particles],       \
    particles->z[0:num_particles],       \
    particles->v_x[0:num_particles],     \
    particles->v_y[0:num_particles],     \
    particles->v_z[0:num_particles],     \
    particles->dp_x[0:num_particles],    \
    particles->dp_y[0:num_particles],    \
    particles->dp_z[0:num_particles],    \
    particles->w_x[0:num_particles],     \
    particles->w_y[0:num_particles],     \
    particles->w_z[0:num_particles],     \
    particles->density[0:num_particles], \
    particles->lambda[0:num_particles],  \
    particles->id[0:num_particles]       \
  )

}

void FinalizeParticles(struct Particles *particles) {
  #pragma acc exit data delete( \
    particles->x_star,          \
    particles->y_star,          \
    particles->z_star,          \
    particles->x,               \
    particles->y,               \
    particles->z,               \
    particles->v_x,             \
    particles->v_y,             \
    particles->v_z,             \
    particles->dp_x,            \
    particles->dp_y,            \
    particles->dp_z,            \
    particles->w_x,             \
    particles->w_y,             \
    particles->w_z,             \
    particles->density,         \
    particles->lambda,          \
    particles->id               \
  )

    free(particles->x_star);
    free(particles->y_star);
    free(particles->z_star);
    free(particles->x);
    free(particles->y);
    free(particles->z);
    free(particles->v_x);
    free(particles->v_y);
    free(particles->v_z);
    free(particles->dp_x);
    free(particles->dp_y);
    free(particles->dp_z);
    free(particles->w_x);
    free(particles->w_y);
    free(particles->w_z);
    free(particles->lambda);
    free(particles->density);
    free(particles->id);
}
