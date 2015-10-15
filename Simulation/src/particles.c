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

  const double W = norm*(h*h - r*r)*(h*h - r*r)*(h*h - r*r);
  return W;
}

#pragma acc routine seq
static inline double DelWPoly(const double r, const double h, const double norm) {
  if(r > h)
    return 0.0;

  const double W = norm*(h*h - r*r)*(h*h - r*r)*r;
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

// http://cg.informatik.uni-freiburg.de/publications/siggraphasia2013/2013_SIGGRAPHASIA_surface_tension_adhesion.pdf
#pragma acc routine seq
static inline double C(const double r, const double h, const double norm) {
  double C;

  if(r > h)
    C = 0.0;
  else if(r <= h/2.0)
    C = (2.0*(h-r)*(h-r)*(h-r)*r*r*r) - (h*h*h*h*h*h/64.0);
  else // h/2 < r < h
    C = (h-r)*(h-r)*(h-r)*r*r*r;

  return norm*C;
}

void ComputeSurfaceTension(struct Particles *restrict particles,
                           const struct Params *restrict params,
                           const struct Neighbors *restrict neighbors) {
  const double rest_density = particles->rest_density;
  const double rest_mass = particles->rest_mass;
  const double h = params->smoothing_radius;
  const double DelW_poly_norm = params->DelW_poly_norm;
  const double dt = params->time_step;
  const int num_particles = particles->local_count;

  // Cohesion Term
  const double C_norm = 32.0/(M_PI*pow(h,9.0));
  const double gama = params->k;

    // OpenACC can't reliably handle SoA...
    const double *restrict x_star = particles->x_star;
    const double *restrict y_star = particles->y_star;
    const double *restrict z_star = particles->z_star;
    double *restrict v_x = particles->v_x;
    double *restrict v_y = particles->v_y;
    double *restrict v_z = particles->v_z;
    double *restrict color_x = particles->color_x;
    double *restrict color_y = particles->color_y;
    double *restrict color_z = particles->color_z;
    const double *restrict density = particles->density;
    const struct NeighborBucket *restrict neighbor_buckets = neighbors->neighbor_buckets;

    //Calculate color for each particle
    #pragma acc parallel loop gang vector vector_length(64)  \
                              present(particles,             \
                              x_star, y_star, z_star,        \
                              color_x, color_y, color_z,     \
                              density,                       \
                              params,                        \
                              neighbors,                     \
                              neighbor_buckets)
    for (int i=0; i<num_particles; ++i) {
      const int p_index = i;
      const struct NeighborBucket *restrict n = &neighbor_buckets[i];

      double color_dx = 0.0;
      double color_dy = 0.0;
      double color_dz = 0.0;

      const double x_star_p = x_star[p_index];
      const double y_star_p = y_star[p_index];
      const double z_star_p = z_star[p_index];

      const double density_p = density[p_index];

      #pragma acc loop seq
      for (int j=0; j<n->count; ++j) {
        const int q_index = n->neighbor_indices[j];
        const double x_diff = x_star_p - x_star[q_index];
        const double y_diff = y_star_p - y_star[q_index];
        const double z_diff = z_star_p - z_star[q_index];

        double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
        if(r_mag < h*0.0001)
          r_mag = h*0.0001;

          const double grad = DelWPoly(r_mag, h, DelW_poly_norm);

          const double mass_q = rest_mass;

          color_dx += grad * mass_q / density[q_index] * x_diff/r_mag;
          color_dy += grad * mass_q / density[q_index] * y_diff/r_mag;
          color_dz += grad * mass_q / density[q_index] * z_diff/r_mag;
      }

      color_x[p_index] = h * color_dx;
      color_y[p_index] = h * color_dy;
      color_z[p_index] = h * color_dz;
    }

    // Calculate surface tension at each particle
    #pragma acc parallel loop gang vector vector_length(64)  \
                              present(particles,             \
                              x_star, y_star, z_star,        \
                              v_x, v_y, v_z,                 \
                              color_x, color_y, color_z,     \
                              density,                       \
                              params,                        \
                              neighbors,                     \
                              neighbor_buckets)
    for (int i=0; i<num_particles; ++i) {
      const int p_index = i;
      const struct NeighborBucket *restrict n = &neighbor_buckets[i];

      double F_surface_x = 0.0;
      double F_surface_y = 0.0;
      double F_surface_z = 0.0;

      const double x_star_p = x_star[p_index];
      const double y_star_p = y_star[p_index];
      const double z_star_p = z_star[p_index];

      const double density_p = density[p_index];

      const double color_x_p = color_x[p_index];
      const double color_y_p = color_y[p_index];
      const double color_z_p = color_z[p_index];

      const double mass_p = rest_mass;

      #pragma acc loop seq
      for (int j=0; j<n->count; ++j) {
        const int q_index = n->neighbor_indices[j];
        const double x_diff = x_star_p - x_star[q_index];
        const double y_diff = y_star_p - y_star[q_index];
        const double z_diff = z_star_p - z_star[q_index];

        double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
        if(r_mag < h*0.0001)
          r_mag = h*0.0001;

         const double mass_q = rest_mass;

          const double c = C(r_mag, h, C_norm);
          const double F_cohesion_x = -20.0*gama * mass_p * mass_q * c * x_diff/r_mag;
          const double F_cohesion_y = -20.0*gama * mass_p * mass_q * c * y_diff/r_mag;
          const double F_cohesion_z = -20.0*gama * mass_p * mass_q * c * z_diff/r_mag;

          const double F_curvature_x = -gama * mass_p * (color_x_p - color_x[q_index]);
          const double F_curvature_y = -gama * mass_p * (color_y_p - color_y[q_index]);
          const double F_curvature_z = -gama * mass_p * (color_z_p - color_z[q_index]);

          const double K = 2.0*rest_density / (density_p + density[q_index]);
          F_surface_x += K * (F_cohesion_x + F_curvature_x);
          F_surface_y += K * (F_cohesion_y + F_curvature_y);
          F_surface_z += K * (F_cohesion_z + F_curvature_z);
      }

      // Mass or density ???
      v_x[p_index] += dt * F_surface_x / density_p;
      v_y[p_index] += dt * F_surface_y / density_p;
      v_z[p_index] += dt * F_surface_z / density_p;
    }
}

void ComputeVorticity(struct Particles *restrict particles,
                      const struct Params *restrict params,
                      const struct Neighbors *restrict neighbors) {

  const double dt = params->time_step;
  const double h = params->smoothing_radius;
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
      if(r_mag < h*0.0001)
        r_mag = h*0.0001;

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
  const double eps = 0.1 * params->smoothing_radius;
  const int num_particles = particles->local_count;
  const double h = params->smoothing_radius;
  const double rest_mass = particles->rest_mass;
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
/*
    const double vort_mag = sqrt(w_x[p_index]*w_x[p_index]
                               + w_y[p_index]*w_y[p_index]
                               + w_z[p_index]*w_z[p_index]);
*/

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
      if(r_mag < h*0.0001)
        r_mag = h*0.0001;

      const double dw = DelW(r_mag, h, DelW_norm);

      const double dw_x = dw*x_diff/r_mag;
      const double dw_y = dw*y_diff/r_mag;
      const double dw_z = dw*z_diff/r_mag;

      // Should this be p_index or q_index???
      const double vort_mag = sqrt(w_x[q_index]*w_x[q_index]
                                 + w_y[q_index]*w_y[q_index]
                                 + w_z[q_index]*w_z[q_index]);

      eta_x += vort_mag*dw_x;
      eta_y += vort_mag*dw_y;
      eta_z += vort_mag*dw_z;
    }

    double eta_mag = sqrt(eta_x*eta_x + eta_y*eta_y + eta_z*eta_z);
    if(eta_mag < h*0.0001)
      eta_mag = h*0.0001;

    const double N_x = eta_x/eta_mag;
    const double N_y = eta_y/eta_mag;
    const double N_z = eta_z/eta_mag;

    const double wx = w_x[p_index];
    const double wy = w_y[p_index];
    const double wz = w_z[p_index];

    const double mass_p = rest_mass;

    v_x[p_index] += dt * eps * (N_y*wz - N_z*wy) / mass_p;
    v_y[p_index] += dt * eps * (N_z*wx - N_x*wz) / mass_p;
    v_z[p_index] += dt * eps * (N_x*wy - N_y*wx) / mass_p;
  }
}

void ApplyViscosity(struct Particles *restrict particles,
                   const struct Params *restrict params,
                   const struct Neighbors *restrict neighbors) {

  const double c = params->c;
  const double h = params->smoothing_radius;
  const double rest_mass = particles->rest_mass;
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
      const double mass_q = rest_mass;
      partial_sum_x += vx_diff * w * mass_q / density[q_index];
      partial_sum_y += vy_diff * w * mass_q / density[q_index];
      partial_sum_z += vz_diff * w * mass_q / density[q_index];
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
  const double rest_mass = particles->rest_mass;
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
    const double mass_p = rest_mass;
    double tmp_density = mass_p * W_0;

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

      const double mass_q = rest_mass;

      const double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
      tmp_density += mass_q * W(r_mag, h, W_norm);
    }

    // Update particle density
    density[p_index] = tmp_density;
  }
}

void ApplyGravity(struct Particles *restrict particles,
                  const struct Params *restrict params) {
  const double dt = params->time_step;
  const double g = -params->g;

  double *v_z = particles->v_z;

  #pragma acc parallel loop present(v_z)
  for (int i=0; i<(particles->local_count); ++i) {
    v_z[i] += g*dt;
  }
}

void UpdatePositionStars(struct Particles *restrict particles,
                         const struct AABB *restrict boundary_global,
                         const struct Obstacle *restrict obstacle,
                         const struct Params *params) {

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
            boundary_global,                   \
            params)
  for (int i=0; i<num_particles; ++i) {
    x_star[i] += dp_x[i];
    y_star[i] += dp_y[i];
    z_star[i] += dp_z[i];

    // Enforce boundary conditions
    ApplyBoundaryConditions(&x_star[i], &y_star[i], &z_star[i], boundary_global,
                            params);
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
  const double rest_density = particles->rest_density;
  const double rest_mass = particles->rest_mass;
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

    double Constraint = density[p_index]/rest_density - 1.0;
    const double Ci = (Constraint < 0.0 ? 0.0 : Constraint);

    double sum_C = 0.0;
    double sum_grad_x = 0.0;
    double sum_grad_y = 0.0;
    double sum_grad_z = 0.0;

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
      if(r_mag < h*0.0001)
        r_mag = h*0.0001;

      const double mass_q = rest_mass;

      const double grad = DelW(r_mag, h, DelW_norm);

      const double grad_x = mass_q * grad*x_diff/r_mag;
      const double grad_y = mass_q * grad*y_diff/r_mag;
      const double grad_z = mass_q * grad*z_diff/r_mag;
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

    // Having particles of different mass would require 1/m_i Term
    // http://mmacklin.com/EG2015PBD.pdf

    sum_C /= (rest_density * rest_density);

    const double epsilon = 0.01; // Jiggle this around orders of magnitude depending on problem
    lambda[i] = -Ci/(sum_C + epsilon);
  }
}

void UpdateDPs(struct Particles *restrict particles,
               const struct Params *restrict params,
               const struct Neighbors *restrict neighbors,
               const int substep) {

  const double h = params->smoothing_radius;
  const double rest_density = particles->rest_density;
  const double rest_mass = particles->rest_mass;
  //const double k = params->k;
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

  const double k_stiff = 0.95;
  const double stiffness = 1.0 - pow((1.0-k_stiff), 1.0/(double)substep);

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
      if(r_mag < h*0.0001)
        r_mag = h*0.0001;

      const double mass_q = rest_mass;
//      const double WdWdq = W(r_mag, h, W_norm)/Wdq;
      const double s_corr = 0.0;//-k * WdWdq*WdWdq*WdWdq*WdWdq;
      const double delW = DelW(r_mag, h, DelW_norm);
      const double dp = mass_q * (lambda[p_index] + lambda[q_index] + s_corr) * delW;

      dpx += dp * x_diff/r_mag;
      dpy += dp * y_diff/r_mag;
      dpz += dp * z_diff/r_mag;
    }

    // Having particles of different mass would require 1/m_i Term
    // http://mmacklin.com/EG2015PBD.pdf

    dp_x[p_index] = stiffness*dpx/rest_density;
    dp_y[p_index] = stiffness*dpy/rest_density;
    dp_z[p_index] = stiffness*dpz/rest_density;
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
                                    boundary_global,        \
                                    params)
  for (int i=0; i<num_particles; ++i) {
    x_star[i] = x[i] + (v_x[i] * dt);
    y_star[i] = y[i] + (v_y[i] * dt);
    z_star[i] = z[i] + (v_z[i] * dt);

    // Enforce boundary conditions before hash
    // Otherwise predicted position can blow up hash
    ApplyBoundaryConditions(&x_star[i], &y_star[i], &z_star[i], boundary_global,
                            params);
//    ApplyBoundaryConditions(&x[i], &y[i], &z[i], boundary_global,
//                            params);
  }

}

// Update particle position and check boundary
void UpdateVelocities(struct Particles *restrict particles,
                      const struct Params *restrict params) {

  const double dt = params->time_step;
  const double v_max = 100.0; // Change to ~sqrt(2*h_max/g) or smoothing_radis/timestep or something

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
                             double *restrict z, const struct AABB *restrict boundary,
                             const struct Params *params) {

  const double h = params->smoothing_radius;
  const double epsion = h*0.0001;

  // Make sure object is not outside boundary
  // The particle must not be equal to boundary max or
  // hash potentially won't pick it up
  // as the particle will in the 'next' after last bin
  if (*x < boundary->min_x)
    *x = boundary->min_x;
  else if(*x  >= boundary->max_x)
    *x = boundary->max_x-epsion;

  if(*y < boundary->min_y)
    *y = boundary->min_y;
  else if(*y  >= boundary->max_y)
    *y = boundary->max_y-epsion;

  if(*z < boundary->min_z)
    *z = boundary->min_z;
  else if(*z >= boundary->max_z)
    *z = boundary->max_z-epsion;
}

void PrintAverageDensity(struct Particles *restrict particles) {

  const int local_count = particles->local_count;
  const double *density = particles->density;

  double average_density = 0.0;
  #pragma acc parallel loop reduction(+:average_density) present(density)
  for (int i=0; i<local_count; ++i) {
    average_density += density[i];
  }

  average_density /= (double)local_count;
  printf("average_density: %.16f\n", average_density);
}

void PrintMaxDensity(struct Particles *restrict particles) {

  const int local_count = particles->local_count;
  const double *density = particles->density;

  double max_density = 0.0;
  #pragma acc parallel loop reduction(max:max_density) present(density)
  for (int i=0; i<local_count; ++i) {
    max_density = fmax(density[i], max_density);
  }

  printf("max density: %.16f\n", max_density);
}

void PrintDensity(struct Particles *restrict particles, const int id) {
  const double *density = particles->density;

  #pragma acc update host(density[id:1])

  printf("density[%d]: %.16f\n", id, density[id]);
}

void PrintAverageDp(struct Particles *restrict particles) {

  const int local_count = particles->local_count;
  const double *dp_x = particles->dp_x;
  const double *dp_y = particles->dp_y;
  const double *dp_z = particles->dp_z;

  double average_dp_x, average_dp_y, average_dp_z = 0.0;
  #pragma acc parallel loop reduction(+:average_dp_x,average_dp_y,average_dp_z) present(dp_x, dp_y, dp_z)
  for (int i=0; i<local_count; ++i) {
    average_dp_x += dp_x[i];
    average_dp_y += dp_y[i];
    average_dp_z += dp_z[i];
  }

  average_dp_x /= (double)local_count;
  average_dp_y /= (double)local_count;
  average_dp_z /= (double)local_count;

  printf("average_dp: (%.16f, %.16f, %.16f)\n", average_dp_x, average_dp_y, average_dp_z);
}

void PrintDp(struct Particles *restrict particles, const int id) {
  const double *dp_x = particles->dp_x;
  const double *dp_y = particles->dp_y;
  const double *dp_z = particles->dp_z;

  #pragma acc update host(dp_x[id:1], dp_y[id:1], dp_z[id:1])

  printf("dp[%d]: (%.16f, %.16f, %.16f)\n", id, dp_x[id], dp_y[id], dp_z[id]);
}

void PrintColor(struct Particles *restrict particles, const int id) {
  const double *color_x = particles->color_x;
  const double *color_y = particles->color_y;
  const double *color_z = particles->color_z;

  #pragma acc update host(color_x[id:1], color_y[id:1], color_z[id:1])

  printf("color[%d]: (%.16f, %.16f, %.16f)\n", id, color_x[id], color_y[id], color_z[id]);
}

void PrintAverageLambda(struct Particles *restrict particles) {

  const int local_count = particles->local_count;
  const double *lambda = particles->lambda;

  double average_lambda = 0.0;
  #pragma acc parallel loop reduction(+:average_lambda) present(lambda)
  for (int i=0; i<local_count; ++i) {
    average_lambda += lambda[i];
  }

  average_lambda /= (double)local_count;

  printf("average_lambda: %.16f\n", average_lambda);
}

void PrintLambda(struct Particles *restrict particles, const int id) {
  const double *lambda = particles->lambda;

  #pragma acc update host(lambda[id:1])

  printf("lambda[%d]: %.16f\n", id, lambda[id]);
}

void PrintAverageV(struct Particles *restrict particles) {

  const int local_count = particles->local_count;
  const double *v_x = particles->v_x;
  const double *v_y = particles->v_y;
  const double *v_z = particles->v_z;

  double average_v_x, average_v_y, average_v_z = 0.0;
  #pragma acc parallel loop reduction(+:average_v_x,average_v_y,average_v_z) present(v_x, v_y, v_z)
  for (int i=0; i<local_count; ++i) {
    average_v_x += v_x[i];
    average_v_y += v_y[i];
    average_v_z += v_z[i];
  }

  average_v_x /= (double)local_count;
  average_v_y /= (double)local_count;
  average_v_z /= (double)local_count;

  printf("average_v: (%.16f, %.16f, %.16f)\n", average_v_x, average_v_y, average_v_z);
}

void PrintVelocity(struct Particles *restrict particles, const int id) {
  const double *v_x = particles->v_x;
  const double *v_y = particles->v_y;
  const double *v_z = particles->v_z;

  #pragma acc update host(v_x[id:1], v_y[id:1], v_z[id:1])

  printf("v[%d]: (%.16f, %.16f, %.16f)\n", id, v_x[id], v_y[id], v_z[id]);
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
  particles->color_x    = SAFE_ALLOC(num_particles, sizeof(double));
  particles->color_y    = SAFE_ALLOC(num_particles, sizeof(double));
  particles->color_z    = SAFE_ALLOC(num_particles, sizeof(double));
  particles->density = SAFE_ALLOC(num_particles, sizeof(double));
  particles->lambda  = SAFE_ALLOC(num_particles, sizeof(double));
  particles->id      = SAFE_ALLOC(num_particles, sizeof(int));

  ConstructFluidVolume(particles, params, fluid_volume_initial);

  particles->rest_density = 1000.0;
  particles->rest_mass = params->fluid_volume * particles->rest_density / (double)particles->global_count;
  // For some reason the initial density estimate is slightly too high
  particles->rest_mass *= .991;

  printf("Rest mass: %f\n", particles->rest_mass);
  printf("Rest Density: %f\n", particles->rest_density);

  // Initialize particle values`
  for (int i=0; i<particles->local_count; ++i) {
    particles->v_x[i] = 0.0;
    particles->v_y[i] = 0.0;
    particles->v_z[i] = 0.0;
    particles->x_star[i] = particles->x[i];
    particles->y_star[i] = particles->y[i];
    particles->z_star[i] = particles->z[i];
    particles->density[i] = particles->rest_density;
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
    particles->color_x[0:num_particles],     \
    particles->color_y[0:num_particles],     \
    particles->color_z[0:num_particles],     \
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
    particles->color_x,             \
    particles->color_y,             \
    particles->color_z,             \
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
  free(particles->color_x);
  free(particles->color_y);
  free(particles->color_z);
  free(particles->lambda);
  free(particles->density);
  free(particles->id);
}
