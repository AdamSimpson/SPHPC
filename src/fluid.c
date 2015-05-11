#include "fluid.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Compiler must support vardiac macros
#ifdef DEBUG
# define DEBUG_PRINT(str, args...) printf("DEBUG: %s:%d:%s(): " str, \
 __FILE__, __LINE__, __func__, ##args)
#else
# define DEBUG_PRINT(str, args...) do {} while (0)
#endif

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

////////////////////////////////////////////////////////////////////////////
// Smoothing Kernels
// Don't use pow as it's POWerfully slow
///////////////////////////////////////////////////////////////////////////

// (h^2 - r^2)^3 normalized in 3D (poly6)
static double W(const double r, const double h)
{
    if(r > h)
        return 0.0;

    const double C = 315.0/(64.0*M_PI* h*h*h*h*h*h*h*h*h);
    const double W = C*(h*h-r*r)*(h*h-r*r)*(h*h-r*r);
    return W;
}

// Gradient (h-r)^3 normalized in 3D (Spikey) magnitude
// Need to multiply by r/|r|
static double del_W(const double r, const double h)
{
    if(r > h)
        return 0.0;

    const double C = -45.0/(M_PI * h*h*h*h*h*h);
    const double del_W = C*(h-r)*(h-r);
    return del_W;
}

////////////////////////////////////////////////////////////////////////////
// Particle attribute computations
////////////////////////////////////////////////////////////////////////////

void vorticity_confinement(FluidParticle *const fluid_particles,
                           const Params *const params,
                           const Neighbors *const neighbors)
{
    const double dt = params->time_step;
    const double eps = 5.0;

    // Calculate vorticy at each particle
    for(int i=0; i<params->number_fluid_particles_local; i++)
    {
        FluidParticle *const p = &fluid_particles[i];
        const Neighbor *const n = &neighbors->particle_neighbors[i];

        double vort_x = 0.0;
        double vort_y = 0.0;
        double vort_z = 0.0;

        for(int j=0; j<n->number_fluid_neighbors; j++)
        {
            const FluidParticle *const q = &fluid_particles[n->neighbor_indices[j]];

            const double x_diff = p->x_star - q->x_star;
            const double y_diff = p->y_star - q->y_star;
            const double z_diff = p->z_star - q->z_star;

            const double vx_diff = q->v_x - p->v_x;
            const double vy_diff = q->v_y - p->v_y;
            const double vz_diff = q->v_z - p->v_z;

            double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            if(r_mag < 0.0001)
              r_mag = 0.0001;

            const double dw = del_W(r_mag, params->smoothing_radius);

            const double dw_x = dw*x_diff/r_mag;
            const double dw_y = dw*y_diff/r_mag;
            const double dw_z = dw*z_diff/r_mag;

            vort_x +=  vy_diff*dw_z - vz_diff*dw_y;
            vort_y +=  vz_diff*dw_x - vx_diff*dw_z;
            vort_z +=  vx_diff*dw_y - vy_diff*dw_x;
        }

        p->w_x = vort_x;
        p->w_y = vort_y;
        p->w_z = vort_z;
    }

    // Apply vorticity confinement
    for(int i=0; i<params->number_fluid_particles_local; i++)
    {
        FluidParticle *const p = &fluid_particles[i];
        const Neighbor *const n = &neighbors->particle_neighbors[i];

        double eta_x  = 0.0;
        double eta_y  = 0.0;
        double eta_z  = 0.0;

        for(int j=0; j<n->number_fluid_neighbors; j++)
        {
            const FluidParticle *const q = &fluid_particles[n->neighbor_indices[j]];

            const double x_diff = p->x_star - q->x_star;
            const double y_diff = p->y_star - q->y_star;
            const double z_diff = p->z_star - q->z_star;

            double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            if(r_mag < 0.0001)
              r_mag = 0.0001;

            const double dw = del_W(r_mag, params->smoothing_radius);

            const double dw_x = dw*x_diff/r_mag;
            const double dw_y = dw*y_diff/r_mag;
            const double dw_z = dw*z_diff/r_mag;

            const double vort_mag = sqrt(q->w_x*q->w_x + q->w_y*q->w_y + q->w_z*q->w_z);

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

        p->v_x += dt*eps * (N_y*p->w_z - N_z*p->w_y);
        p->v_y += dt*eps * (N_z*p->w_x - N_x*p->w_z);
        p->v_z += dt*eps * (N_x*p->w_y - N_y*p->w_x);

    }
}

void XSPH_viscosity(FluidParticle *const fluid_particles,
                    const Params *const params,
                    const Neighbors *const neighbors)
{
    const double c = params->c;
    const double h = params->smoothing_radius;

    for(int i=0; i<params->number_fluid_particles_local; i++)
    {
        const Neighbor *const n = &neighbors->particle_neighbors[i];

        double partial_sum_x = 0.0;
        double partial_sum_y = 0.0;
        double partial_sum_z = 0.0;

        for(int j=0; j<n->number_fluid_neighbors; j++)
        {
            const int q_index = n->neighbor_indices[j];
            const double x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            const double y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            const double z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            const double vx_diff = fluid_particles[q_index].v_x - fluid_particles[i].v_x;
            const double vy_diff = fluid_particles[q_index].v_y - fluid_particles[i].v_y;
            const double vz_diff = fluid_particles[q_index].v_z - fluid_particles[i].v_z;

            const double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff * z_diff*z_diff);

            const double w = W(r_mag, h);

            // http://mmacklin.com/pbf_sig_preprint.pdf is missing 1/sigma contribution
            // see: http://www.cs.ubc.ca/~rbridson/docs/schechter-siggraph2012-ghostsph.pdf
            partial_sum_x += vx_diff * w / fluid_particles[q_index].density;
            partial_sum_y += vy_diff * w / fluid_particles[q_index].density;
            partial_sum_z += vz_diff * w / fluid_particles[q_index].density;
        }

        partial_sum_x *= c;
        partial_sum_y *= c;
        partial_sum_z *= c;

        fluid_particles[i].v_x += partial_sum_x;
        fluid_particles[i].v_y += partial_sum_y;
        fluid_particles[i].v_z += partial_sum_z;
    }
}

void compute_densities(FluidParticle *const fluid_particles,
                       const Params *const params,
                       const Neighbors *const neighbors)
{
    const double h = params->smoothing_radius;

    for(int i=0; i<params->number_fluid_particles_local; i++)
    {
        const Neighbor *const n = &neighbors->particle_neighbors[i];

        // Own contribution to density
        double density = W(0.0, h);

        // Neighbor contribution
        for(int j=0; j<n->number_fluid_neighbors; j++)
        {
            const int q_index = n->neighbor_indices[j];
            const double x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            const double y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            const double z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            const double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            density += W(r_mag, h);
        }

        // Update particle density
        fluid_particles[i].density = density;
    }
}

void apply_gravity(FluidParticle *const fluid_particles,
                   const Params *const params)
{
    const double dt = params->time_step;
    const double g = -params->g;

    for(int i=0; i<(params->number_fluid_particles_local); i++) {
        fluid_particles[i].v_y += g*dt;
    }
}

void update_dp_positions(FluidParticle *const fluid_particles,
                         const Params *const params,
                         const AABB *const boundary_global)
{
    for(int i=0; i<(params->number_fluid_particles_local); i++) {

        fluid_particles[i].x_star += fluid_particles[i].dp_x;
        fluid_particles[i].y_star += fluid_particles[i].dp_y;
        fluid_particles[i].z_star += fluid_particles[i].dp_z;

        // Enforce boundary conditions
        boundary_conditions(fluid_particles, i, boundary_global);

    }
}

void update_positions(FluidParticle *const fluid_particles,
                      const Params *const params)
{
     for(int i=0; i<(params->number_fluid_particles_local); i++) {
        fluid_particles[i].x = fluid_particles[i].x_star;
        fluid_particles[i].y = fluid_particles[i].y_star;
        fluid_particles[i].z = fluid_particles[i].z_star;
    }
}

void calculate_lambda(FluidParticle *const fluid_particles,
                      const Params *const params,
                      const Neighbors *const neighbors)
{
    for(int i=0; i<params->number_fluid_particles_local; i++)
    {
        const Neighbor *const n = &neighbors->particle_neighbors[i];

        const double Ci = fluid_particles[i].density/params->rest_density - 1.0;

        double sum_C = 0.0;
        double sum_grad_x = 0.0;
        double sum_grad_y = 0.0;
        double sum_grad_z = 0.0;

        for(int j=0; j<n->number_fluid_neighbors; j++)
        {
            const int q_index = n->neighbor_indices[j];
            const double x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            const double y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            const double z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            if(r_mag < 0.0001) {
              r_mag = 0.0001;
              DEBUG_PRINT("pstar: %f, %f, %f grad: %f\n", fluid_particles[i].x_star, fluid_particles[i].y_star,fluid_particles[i].z_star, grad);
            }

            const double grad = del_W(r_mag, params->smoothing_radius);

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
        fluid_particles[i].lambda = -Ci/(sum_C + epsilon);
    }
}

void update_dp(FluidParticle *const fluid_particles,
               const Params *const params,
               const Neighbors *const neighbors)
{
    const double h = params->smoothing_radius;
    const double k = params->k;
    const double dq = params->dq;
    const double Wdq = W(dq, h);

    for(int i=0; i<params->number_fluid_particles_local; i++)
    {
        const Neighbor *const n = &neighbors->particle_neighbors[i];

        double dp_x = 0.0;
        double dp_y = 0.0;
        double dp_z = 0.0;

        for(int j=0; j<n->number_fluid_neighbors; j++)
        {
            const int q_index = n->neighbor_indices[j];
            const double x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            const double y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            const double z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            double r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            if(r_mag < 0.0001)
              r_mag = 0.0001;

            const double WdWdq = W(r_mag, h)/Wdq;
            const double s_corr = -k*WdWdq*WdWdq*WdWdq*WdWdq;
            const double dp = (fluid_particles[i].lambda + fluid_particles[q_index].lambda + s_corr)*del_W(r_mag, h);

            dp_x += dp*x_diff/r_mag;
            dp_y += dp*y_diff/r_mag;
            dp_z += dp*z_diff/r_mag;
        }
        fluid_particles[i].dp_x = dp_x/params->rest_density;
        fluid_particles[i].dp_y = dp_y/params->rest_density;
        fluid_particles[i].dp_z = dp_z/params->rest_density;
    }
}

// Identify out of bounds particles and send them to appropriate rank
void identify_oob_particles(FluidParticle *const fluid_particles,
                            Params *const params,
                            Communication *const communication)
{
    OOB *const out_of_bounds = &communication->out_of_bounds;

    // Reset OOB numbers
    out_of_bounds->number_oob_particles_left = 0;
    out_of_bounds->number_oob_particles_right = 0;

    for(int i=0; i<params->number_fluid_particles_local; i++) {
        // Set OOB particle indices and update number
        if (fluid_particles[i].x_star < params->node_start_x)
            out_of_bounds->oob_indices_left[out_of_bounds->number_oob_particles_left++] = i;
        else if (fluid_particles[i].x_star > params->node_end_x)
            out_of_bounds->oob_indices_right[out_of_bounds->number_oob_particles_right++] = i;
    }

   // Transfer particles that have left the processor bounds
   transfer_OOB_particles(communication, fluid_particles, params);
}

// Predict position
void predict_positions(FluidParticle *const fluid_particles,
                       const Params *const params,
                       const AABB *const boundary_global)
{
    const double dt = params->time_step;

    for(int i=0; i<params->number_fluid_particles_local; i++) {
        fluid_particles[i].x_star = fluid_particles[i].x + (fluid_particles[i].v_x * dt);
        fluid_particles[i].y_star = fluid_particles[i].y + (fluid_particles[i].v_y * dt);
        fluid_particles[i].z_star = fluid_particles[i].z + (fluid_particles[i].v_z * dt);

        // Enforce boundary conditions before hash
        // Otherwise predicted position can blow up hash
        boundary_conditions(fluid_particles, i, boundary_global);
    }
}

void check_velocity(double *const v_x, double *const v_y, double *const v_z)
{
    const double v_max = 30.0;

    if(*v_x > v_max)
        *v_x = v_max;
    else if(*v_x < -v_max)
        *v_x = -v_max;
    if(*v_y > v_max)
        *v_y = v_max;
    else if(*v_y < -v_max)
        *v_y = -v_max;
    if(*v_z > v_max)
        *v_z = v_max;
    else if(*v_z < -v_max)
        *v_z = -v_max;

}

// Update particle position and check boundary
void update_velocities(FluidParticle *const fluid_particles,
                       const Params *const params)
{
    const double dt = params->time_step;

    // Update local and halo particles, update halo so that XSPH visc. is correct
    for(int i=0; i<params->number_fluid_particles_local + params->number_halo_particles_left + params->number_halo_particles_right; i++) {
        double v_x = (fluid_particles[i].x_star - fluid_particles[i].x)/dt;
        double v_y = (fluid_particles[i].y_star - fluid_particles[i].y)/dt;
        double v_z = (fluid_particles[i].z_star - fluid_particles[i].z)/dt;

        check_velocity(&v_x, &v_y, &v_z);

        fluid_particles[i].v_x = v_x;
        fluid_particles[i].v_y = v_y;
        fluid_particles[i].v_z = v_z;
    }
}

// Assume AABB with min point being axis origin
void boundary_conditions(FluidParticle *const fluid_particles,
                         const unsigned int i,
                         const AABB *const boundary)
{

    // Make sure object is not outside boundary
    // The particle must not be equal to boundary max or hash potentially won't pick it up
    // as the particle will in the 'next' after last bin
    if(fluid_particles[i].x_star  < boundary->min_x)
        fluid_particles[i].x_star = boundary->min_x;
    else if(fluid_particles[i].x_star  > boundary->max_x)
        fluid_particles[i].x_star = boundary->max_x-0.00001;

    if(fluid_particles[i].y_star  <  boundary->min_y)
        fluid_particles[i].y_star = boundary->min_y;
    else if(fluid_particles[i].y_star  > boundary->max_y)
        fluid_particles[i].y_star = boundary->max_y-0.00001;

    if(fluid_particles[i].z_star  <  boundary->min_z)
        fluid_particles[i].z_star = boundary->min_z;
    else if(fluid_particles[i].z_star  > boundary->max_z)
        fluid_particles[i].z_star = boundary->max_z-0.00001;
}

void init_particles(FluidParticle *const fluid_particles,
                    Params *const params,
                    const AABB *const water)
{
    // Create fluid volume
    construct_fluid_volume(fluid_particles, params, water);

    // Initialize particle values
    for(int i=0; i<params->number_fluid_particles_local; i++) {
        fluid_particles[i].x_star = fluid_particles[i].x;
        fluid_particles[i].y_star = fluid_particles[i].y;
        fluid_particles[i].z_star = fluid_particles[i].z;
        fluid_particles[i].v_x = 0.0;
        fluid_particles[i].v_y = 0.0;
        fluid_particles[i].v_z = 0.0;
        fluid_particles[i].dp_x = 0.0;
        fluid_particles[i].dp_y = 0.0;
        fluid_particles[i].dp_z = 0.0;
        fluid_particles[i].w_x = 0.0;
        fluid_particles[i].w_y = 0.0;
        fluid_particles[i].w_z = 0.0;
        fluid_particles[i].lambda = 0.0;
        fluid_particles[i].density = params->rest_density;
    }
}

void allocate_fluid(FluidParticle **fluid_particles,
                    const Params *const params)
{
    *fluid_particles = calloc(params->max_fluid_particles_local, sizeof(FluidParticle));
    if(*fluid_particles == NULL)
        printf("Could not allocate fluid_particles\n");
}
