#include "particles.h"
#include "simulation.h"
#include "communication.h"
#include "neighbors.h"
#include "geometry.h"
#include "obstacle.h"
#include "print_macros.h"
#include "safe_alloc.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "openacc.h"

///////////////////////////////////////////////////////////////////////////
// Smoothing Kernels
///////////////////////////////////////////////////////////////////////////

// (h^2 - r^2)^3 normalized in 3D (poly6)
#pragma acc routine seq
static inline double W(const double r, const double h, const double norm) {
  if(r > h)
    return 0.0;

  const double W = norm*(h*h-r*r)*(h*h-r*r)*(h*h-r*r);
  return W;
}

// Gradient (h-r)^3 normalized in 3D (Spikey) magnitude
// Need to multiply by r/|r|
#pragma acc routine seq
static inline double DelW(const double r, const double h, const double norm) {
  if(r > h)
    return 0.0;

  const double DelW = norm*(h-r)*(h-r);
  return DelW;
}

////////////////////////////////////////////////////////////////////////////
// Particle attribute computations
////////////////////////////////////////////////////////////////////////////

void ComputeVorticity(struct Particles *restrict particles,
                      const struct Params *restrict params,
                      const struct Neighbors *restrict neighbors) {

  const double dt = params->time_step;
  const double h = params->smoothing_radius;
  const double eps = 0.001 * params->smoothing_radius;
  const double DelW_norm = params->DelW_norm;
  const int num_particles = particles->local_count;

  // OpenACC can't reliably handle SoA...
  const double *restrict x_star = particles->x_star;
  const double *restrict y_star = particles->y_star;
  const double *restrict z_star = particles->z_star;
  const double *restrict v_x = particles->v_x;
  const double *restrict v_y = particles->v_y;
  const double *restrict v_z = particles->v_z;
  double *restrict w_x = particles->w_x;
  double *restrict w_y = particles->w_y;
  double *restrict w_z = particles->w_z;
  const struct NeighborBucket *restrict neighbor_buckets = neighbors->neighbor_buckets;

  // Calculate vorticy at each particle
  #pragma acc parallel loop gang vector vector_length(64)  \
                            present(particles,             \
                            x_star, y_star, z_star,        \
                            v_x, v_y, v_z,                 \
                            w_x, w_y, w_z,                 \
                            params,                        \
                            neighbors,                     \
                            neighbor_buckets)
  for (int i=0; i<num_particles; ++i) {
    const int p_index = i;
    const struct NeighborBucket *restrict n = &neighbor_buckets[i];

    double vort_x = 0.0;
    double vort_y = 0.0;
    double vort_z = 0.0;

    const double x_star_p = x_star[p_index];
    const double y_star_p = y_star[p_index];
    const double z_star_p = z_star[p_index];

    const double v_x_p = v_x[p_index];
    const double v_y_p = v_y[p_index];
    const double v_z_p = v_z[p_index];

    #pragma acc loop seq
    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];
      const double x_diff = x_star_p - x_star[q_index];
      const double y_diff = y_star_p - y_star[q_index];
      const double z_diff = z_star_p - z_star[q_index];

      const double vx_diff = v_x[q_index] - v_x_p;
      const double vy_diff = v_y[q_index] - v_y_p;
      const double vz_diff = v_z[q_index] - v_z_p;

      double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
      if(r_mag < 0.0001)
        r_mag = 0.0001;

        const double dw = DelW(r_mag, h, DelW_norm);

        const double dw_x = dw*x_diff/r_mag;
        const double dw_y = dw*y_diff/r_mag;
        const double dw_z = dw*z_diff/r_mag;

        vort_x +=  vy_diff*dw_z - vz_diff*dw_y;
        vort_y +=  vz_diff*dw_x - vx_diff*dw_z;
        vort_z +=  vx_diff*dw_y - vy_diff*dw_x;
    }

    w_x[p_index] = vort_x;
    w_y[p_index] = vort_y;
    w_z[p_index] = vort_z;
  }
}

void ApplyVorticityConfinement(struct Particles *restrict particles,
                               const struct Params *restrict params,
                               const struct Neighbors *restrict neighbors) {
  const double dt = params->time_step;
  const double eps = 0.001 * params->smoothing_radius;
  const int num_particles = particles->local_count;
  const double h = params->smoothing_radius;
  const double DelW_norm = params->DelW_norm;

  // OpenACC can't reliably handle SoA...
  const double *restrict x_star = particles->x_star;
  const double *restrict y_star = particles->y_star;
  const double *restrict z_star = particles->z_star;
  double *restrict v_x = particles->v_x;
  double *restrict v_y = particles->v_y;
  double *restrict v_z = particles->v_z;
  const double *restrict w_x = particles->w_x;
  const double *restrict w_y = particles->w_y;
  const double *restrict w_z = particles->w_z;
  const struct NeighborBucket *restrict neighbor_buckets = neighbors->neighbor_buckets;

  // Apply vorticity confinement
  #pragma acc parallel loop gang vector vector_length(64) \
                            present(particles,            \
                            x_star, y_star, z_star,       \
                            v_x, v_y, v_z,                \
                            w_x, w_y, w_z,                \
                            params,                       \
                            neighbors,                    \
                            neighbor_buckets)
  for (int i=0; i<particles->local_count; ++i) {
    const int p_index = i;
    const struct NeighborBucket *restrict n = &neighbor_buckets[i];

    double eta_x  = 0.0;
    double eta_y  = 0.0;
    double eta_z  = 0.0;

    // Should this be p_index or q_index???
    const double vort_mag = sqrt(w_x[p_index]*w_x[p_index]
                               + w_y[p_index]*w_y[p_index]
                               + w_z[p_index]*w_z[p_index]);

    const double x_star_p = x_star[p_index];
    const double y_star_p = y_star[p_index];
    const double z_star_p = z_star[p_index];

    #pragma acc loop seq
    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];

      const double x_diff = x_star_p - x_star[q_index];
      const double y_diff = y_star_p - y_star[q_index];
      const double z_diff = z_star_p - z_star[q_index];

      double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
      if(r_mag < 0.0001)
        r_mag = 0.0001;

      const double dw = DelW(r_mag, h, DelW_norm);

      const double dw_x = dw*x_diff/r_mag;
      const double dw_y = dw*y_diff/r_mag;
      const double dw_z = dw*z_diff/r_mag;

      // Should this be p_index or q_index???
      /*
      const double vort_mag = sqrt(w_x[q_index]*w_x[q_index]
                                 + w_y[q_index]*w_y[q_index]
                                 + w_z[q_index]*w_z[q_index]);
      */

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

    const double wx = w_x[p_index];
    const double wy = w_y[p_index];
    const double wz = w_z[p_index];

    v_x[p_index] += eps * (N_y*wz - N_z*wy);
    v_y[p_index] += eps * (N_z*wx - N_x*wz);
    v_z[p_index] += eps * (N_x*wy - N_y*wx);
  }
}

void ApplyViscosity(struct Particles *restrict particles,
                   const struct Params *restrict params,
                   const struct Neighbors *restrict neighbors) {

  const double c = params->c;
  const double h = params->smoothing_radius;
  const double W_norm = params->W_norm;

  const int num_particles = particles->local_count;

  // OpenACC can't reliably handle SoA...
  const double *restrict x_star = particles->x_star;
  const double *restrict y_star = particles->y_star;
  const double *restrict z_star = particles->z_star;
  double *restrict v_x = particles->v_x;
  double *restrict v_y = particles->v_y;
  double *restrict v_z = particles->v_z;
  const double *restrict density = particles->density;
  const struct NeighborBucket *restrict neighbor_buckets = neighbors->neighbor_buckets;

  #pragma acc parallel loop gang vector vector_length(64) \
                            present(particles,            \
                            x_star, y_star, z_star,       \
                            v_x, v_y, v_z,                \
                            density,                      \
                            params,                       \
                            neighbors,                    \
                            neighbor_buckets)
  for (int i=0; i < num_particles; ++i) {
    const int p_index = i;
    const struct NeighborBucket *restrict n = &neighbor_buckets[i];

    double partial_sum_x = 0.0;
    double partial_sum_y = 0.0;
    double partial_sum_z = 0.0;

    const double x_star_p = x_star[p_index];
    const double y_star_p = y_star[p_index];
    const double z_star_p = z_star[p_index];

    const double v_x_p = v_x[p_index];
    const double v_y_p = v_y[p_index];
    const double v_z_p = v_z[p_index];

    #pragma acc loop seq
    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];
      const double x_diff = x_star_p - x_star[q_index];
      const double y_diff = y_star_p - y_star[q_index];
      const double z_diff = z_star_p - z_star[q_index];

      const double vx_diff = v_x[q_index] - v_x_p;
      const double vy_diff = v_y[q_index] - v_y_p;
      const double vz_diff = v_z[q_index] - v_z_p;

      const double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff * z_diff*z_diff);

      const double w = W(r_mag, h, W_norm);

      // http://mmacklin.com/pbf_sig_preprint.pdf is missing 1/sigma contribution
      // see: http://www.cs.ubc.ca/~rbridson/docs/schechter-siggraph2012-ghostsph.pdf
      partial_sum_x += vx_diff * w / density[q_index];
      partial_sum_y += vy_diff * w / density[q_index];
      partial_sum_z += vz_diff * w / density[q_index];
    }

    partial_sum_x *= c;
    partial_sum_y *= c;
    partial_sum_z *= c;

    v_x[p_index] += partial_sum_x;
    v_y[p_index] += partial_sum_y;
    v_z[p_index] += partial_sum_z;
  }

}

void ComputeDensities(struct Particles *restrict particles,
                      const struct Params *restrict params,
                      const struct Neighbors *restrict neighbors) {

  const double h = params->smoothing_radius;
  const double W_norm = params->W_norm;
  const double W_0 = W(0.0, h, W_norm);
  const int num_particles = particles->local_count;

  // OpenACC can't reliably handle SoA...
  const double *restrict x_star = particles->x_star;
  const double *restrict y_star = particles->y_star;
  const double *restrict z_star = particles->z_star;
  double *restrict density = particles->density;
  const struct NeighborBucket *restrict neighbor_buckets = neighbors->neighbor_buckets;

  #pragma acc parallel loop gang vector vector_length(64) \
                                        present(particles, \
                                        x_star, y_star, z_star,            \
                                        density,                           \
                                        params,                            \
                                        neighbors,                         \
                                        neighbor_buckets)
    for (int i=0; i<particles->local_count; ++i) {
    const int p_index = i;
    const struct NeighborBucket *restrict n = &neighbor_buckets[i];

    // Own contribution to density
    double tmp_density = W_0;

    const double x_star_p = x_star[p_index];
    const double y_star_p = y_star[p_index];
    const double z_star_p = z_star[p_index];

    // Neighbor contribution
    #pragma acc loop seq
    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];
      const double x_diff = x_star_p - x_star[q_index];
      const double y_diff = y_star_p - y_star[q_index];
      const double z_diff = z_star_p - z_star[q_index];

      const double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
      tmp_density += W(r_mag, h, W_norm);
    }

    // Update particle density
    density[p_index] = tmp_density;
  }
}

void ApplyGravity(struct Particles *restrict particles,
                  const struct Params *restrict params) {
  const double dt = params->time_step;
  const double g = -params->g;

  double *v_y = particles->v_y;

  #pragma acc parallel loop present(v_y)
  for (int i=0; i<(particles->local_count); ++i) {
    v_y[i] += g*dt;
  }
}

void UpdatePositionStars(struct Particles *restrict particles,
                         const struct AABB *restrict boundary_global,
                         const struct Obstacle *restrict obstacle) {

  const int num_particles = particles->local_count;

  // OpenACC can't reliably handle SoA...
  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;
  double *dp_x = particles->dp_x;
  double *dp_y = particles->dp_y;
  double *dp_z = particles->dp_z;

  #pragma acc parallel loop present(particles, \
            x_star, y_star, z_star,            \
            dp_x, dp_y, dp_z,                  \
            boundary_global)
  for (int i=0; i<num_particles; ++i) {
    x_star[i] += dp_x[i];
    y_star[i] += dp_y[i];
    z_star[i] += dp_z[i];

    // Enforce boundary conditions
    ApplyBoundaryConditions(&x_star[i], &y_star[i], &z_star[i], boundary_global);
  }

  #pragma acc host_data use_device(x_star, y_star, z_star)
  {
    CalculateParticleCollisions(x_star, y_star, z_star, particles->local_count, *obstacle);
  }

}

void UpdatePositions(struct Particles *restrict particles) {

  const int num_particles = particles->local_count;

  // OpenACC can't reliably handle SoA...
  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;
  double *x = particles->x;
  double *y = particles->y;
  double *z = particles->z;

  #pragma acc parallel loop present(particles,  \
                                    x_star, y_star, z_star,  \
                                    x, y, z)
  for (int i=0; i<num_particles; ++i) {
    x[i] = x_star[i];
    y[i] = y_star[i];
    z[i] = z_star[i];
  }
}

void ComputeLambda(struct Particles *restrict particles,
                   const struct Params *restrict params,
                   const struct Neighbors *restrict neighbors) {

  const int num_particles = particles->local_count;
  const double rest_density = params->rest_density;
  const double h = params->smoothing_radius;
  const double DelW_norm = params->DelW_norm;

  // OpenACC can't reliably handle SoA...
  const double *restrict x_star = particles->x_star;
  const double *restrict y_star = particles->y_star;
  const double *restrict z_star = particles->z_star;
  const double *restrict density = particles->density;
  double *restrict lambda = particles->lambda;
  const struct NeighborBucket *restrict neighbor_buckets = neighbors->neighbor_buckets;

  #pragma acc parallel loop gang vector vector_length(64)           \
                                    present(x_star, y_star, z_star, \
                                    density, lambda,                \
                                    params,                         \
                                    neighbors,                      \
                                    neighbor_buckets )
  for (int i=0; i<num_particles; ++i) {
    const int p_index = i;
    const struct NeighborBucket *restrict n = &neighbor_buckets[i];
    const double Ci = density[p_index]/rest_density - 1.0;
    double sum_C = 0.0;
    double sum_grad_x = 0.0;
    double sum_grad_y = 0.0;
    double sum_grad_z = 0.0;

    const double x_star_p = x_star[p_index];
    const double y_star_p = y_star[p_index];
    const double z_star_p = z_star[p_index];

    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];
      const double x_diff = x_star_p - x_star[q_index];
      const double y_diff = y_star_p - y_star[q_index];
      const double z_diff = z_star_p - z_star[q_index];

      double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
      if (r_mag < 0.0001)
        r_mag = 0.0001;

      const double grad = DelW(r_mag, h, DelW_norm);

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
    sum_C += (sum_grad_x*sum_grad_x
            + sum_grad_y*sum_grad_y
            + sum_grad_z*sum_grad_z);

    sum_C *= (1.0/(rest_density * rest_density));

    const double epsilon = 1.0;
    lambda[i] = -Ci/(sum_C + epsilon);
  }
}

void UpdateDPs(struct Particles *restrict particles,
               const struct Params *restrict params,
               const struct Neighbors *restrict neighbors) {

  const double h = params->smoothing_radius;
  const double rest_density = params->rest_density;
  const double k = params->k;
  const double W_norm = params->W_norm;
  const double DelW_norm = params->DelW_norm;
  const double Wdq = W(params->dq, h, W_norm);
  const int num_particles = particles->local_count;

  // OpenACC can't reliably handle SoA...
  // PGI should use read only data cache if const restrict is used
  const double *restrict x_star  = particles->x_star;
  const double *restrict y_star  = particles->y_star;
  const double *restrict z_star  = particles->z_star;
  double *restrict dp_x    = particles->dp_x;
  double *restrict dp_y    = particles->dp_y;
  double *restrict dp_z    = particles->dp_z;
  const double *restrict lambda = particles->lambda;
  const struct NeighborBucket *restrict neighbor_buckets = neighbors->neighbor_buckets;

  // acc kernels seems to work just as well
  #pragma acc parallel loop gang vector vector_length(64)   \
                                    present(particles,      \
                                    x_star, y_star, z_star, \
                                    dp_x, dp_y, dp_z,       \
                                    lambda,                 \
                                    params,                 \
                                    neighbors,              \
                                    neighbor_buckets )
  for (int i=0; i<num_particles; ++i) {
    const int p_index = i;
    const struct NeighborBucket *restrict n = &neighbor_buckets[i];

    double dpx = 0.0;
    double dpy = 0.0;
    double dpz = 0.0;

    const double x_star_p = x_star[p_index];
    const double y_star_p = y_star[p_index];
    const double z_star_p = z_star[p_index];

    #pragma acc loop seq
    for (int j=0; j<n->count; ++j) {
      const int q_index = n->neighbor_indices[j];
      const double x_diff = x_star_p - x_star[q_index];
      const double y_diff = y_star_p - y_star[q_index];
      const double z_diff = z_star_p - z_star[q_index];

      double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
      if(r_mag < 0.0001)
        r_mag = 0.0001;

      const double WdWdq = W(r_mag, h, W_norm)/Wdq;
      const double delW = DelW(r_mag, h, DelW_norm);
      const double s_corr = -k * WdWdq*WdWdq*WdWdq*WdWdq;
      const double dp = (lambda[p_index] + lambda[q_index] + s_corr) * delW;

      dpx += dp * x_diff/r_mag;
      dpy += dp * y_diff/r_mag;
      dpz += dp * z_diff/r_mag;
    }

    dp_x[p_index] = dpx/rest_density;
    dp_y[p_index] = dpy/rest_density;
    dp_z[p_index] = dpz/rest_density;
  }
}

void PredictPositions(struct Particles *restrict particles,
                      const struct Params *restrict params,
                      const struct AABB *restrict boundary_global) {
  const double dt = params->time_step;

  const int num_particles = particles->local_count;

  // OpenACC can't reliably handle SoA...
  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;
  double *x      = particles->x;
  double *y      = particles->y;
  double *z      = particles->z;
  double *v_x    = particles->v_x;
  double *v_y    = particles->v_y;
  double *v_z    = particles->v_z;

  #pragma acc parallel loop present(x_star, y_star, z_star, \
                                    x, y, z,                \
                                    v_x, v_y, v_z,          \
                                    boundary_global)
  for (int i=0; i<num_particles; ++i) {
    x_star[i] = x[i] + (v_x[i] * dt);
    y_star[i] = y[i] + (v_y[i] * dt);
    z_star[i] = z[i] + (v_z[i] * dt);

    // Enforce boundary conditions before hash
    // Otherwise predicted position can blow up hash
    ApplyBoundaryConditions(&x_star[i], &y_star[i], &z_star[i], boundary_global);
  }

}

// Update particle position and check boundary
void UpdateVelocities(struct Particles *restrict particles,
                      const struct Params *restrict params) {

  const double dt = params->time_step;
  const double v_max = 30.0;

  // Update local and halo particles, update halo so that XSPH visc. is correct
  // NEED TO RETHINK THIS
  int total_particles = particles->local_count
                      + particles->halo_count_left
                      + particles->halo_count_right;

  // OpenACC can't reliably handle SoA...
  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;
  double *x      = particles->x;
  double *y      = particles->y;
  double *z      = particles->z;
  double *v_x    = particles->v_x;
  double *v_y    = particles->v_y;
  double *v_z    = particles->v_z;

  #pragma acc parallel loop present(particles,              \
                                    x_star, y_star, z_star, \
                                    x, y, z,                \
                                    v_x, v_y, v_z)
  for(int i=0; i<total_particles; ++i) {
    double vx = (x_star[i] - x[i])/dt;
    double vy = (y_star[i] - y[i])/dt;
    double vz = (z_star[i] - z[i])/dt;

   if (vx > v_max)
      vx = v_max;
   else if (vx < -v_max)
      vx = -v_max;

   if (vy > v_max)
      vy = v_max;
   else if (vy < -v_max)
      vy = -v_max;

   if (vz > v_max)
      vz = v_max;
   else if (vz < -v_max)
      vz = -v_max;

    v_x[i] = vx;
    v_y[i] = vy;
    v_z[i] = vz;
  }
}

// Assume AABB with min point being axis origin
#pragma acc routine seq
void ApplyBoundaryConditions(double *restrict x, double *restrict y,
                             double *restrict z, const struct AABB *restrict boundary) {

  // Make sure object is not outside boundary
  // The particle must not be equal to boundary max or
  // hash potentially won't pick it up
  // as the particle will in the 'next' after last bin
  if (*x < boundary->min_x)
    *x = boundary->min_x;
  else if(*x  > boundary->max_x)
    *x = boundary->max_x-0.00001;

  if(*y < boundary->min_y)
    *y = boundary->min_y;
  else if(*y  > boundary->max_y)
    *y = boundary->max_y-0.00001;

  if(*z < boundary->min_z)
    *z = boundary->min_z;
  else if(*z > boundary->max_z)
    *z = boundary->max_z-0.00001;
}

void PrintAverageDensity(struct Particles *restrict particles) {
  double average_density = 0.0;

  for (int i=0; i<particles->local_count; ++i) {
    average_density += particles->density[i];
  }
  average_density /= (double)particles->local_count;
  printf("average_density: %.16f\n", average_density);
}

void AllocInitParticles(struct Particles *restrict particles,
                        struct Params *restrict params,
                        struct AABB *restrict fluid_volume_initial) {

  const size_t num_particles = particles->max_local;

  particles->x_star = SAFE_ALLOC(num_particles, sizeof(double));
  particles->y_star = SAFE_ALLOC(num_particles, sizeof(double));
  particles->z_star = SAFE_ALLOC(num_particles, sizeof(double));
  particles->x      = SAFE_ALLOC(num_particles, sizeof(double));
  particles->y      = SAFE_ALLOC(num_particles, sizeof(double));
  particles->z      = SAFE_ALLOC(num_particles, sizeof(double));
  particles->v_x    = SAFE_ALLOC(num_particles, sizeof(double));
  particles->v_y    = SAFE_ALLOC(num_particles, sizeof(double));
  particles->v_z    = SAFE_ALLOC(num_particles, sizeof(double));
  particles->dp_x   = SAFE_ALLOC(num_particles, sizeof(double));
  particles->dp_y   = SAFE_ALLOC(num_particles, sizeof(double));
  particles->dp_z   = SAFE_ALLOC(num_particles, sizeof(double));
  particles->w_x    = SAFE_ALLOC(num_particles, sizeof(double));
  particles->w_y    = SAFE_ALLOC(num_particles, sizeof(double));
  particles->w_z    = SAFE_ALLOC(num_particles, sizeof(double));
  particles->density = SAFE_ALLOC(num_particles, sizeof(double));
  particles->lambda  = SAFE_ALLOC(num_particles, sizeof(double));
  particles->id      = SAFE_ALLOC(num_particles, sizeof(int));

  ConstructFluidVolume(particles, params, fluid_volume_initial);

  // Initialize particle values
  for (int i=0; i<particles->local_count; ++i) {
    particles->v_x[i] = 0.0;
    particles->v_y[i] = 0.0;
    particles->v_z[i] = 0.0;
    particles->x_star[i] = particles->x[i];
    particles->y_star[i] = particles->y[i];
    particles->z_star[i] = particles->z[i];
    particles->density[i] = params->rest_density;
  }

  particles->halo_count_left = 0;
  particles->halo_count_right = 0;

  #pragma acc enter data copyin(particles[:1], \
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

void FinalizeParticles(struct Particles *restrict particles) {
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
