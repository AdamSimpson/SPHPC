#ifndef fluid_fluid_h
#define fluid_fluid_h

typedef struct FLUID_PARTICLE fluid_particle_t;

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
double W(double r, double h);
double del_W(double r, double h);
void XSPH_viscosity(fluid_particle_t *fluid_particles, neighbors_t* neighbors, param_t *params);
void vorticity_confinement(fluid_particle_t *fluid_particles, neighbors_t* neighbors, param_t *params);
void compute_densities(fluid_particle_t *fluid_particles, neighbors_t *neighbors, param_t *params);
void apply_gravity(fluid_particle_t *fluid_particles, param_t *params);
void update_dp_positions(fluid_particle_t *fluid_particles, AABB_t *boundary_global, param_t *params);
void update_positions(fluid_particle_t *fluid_particles, param_t *params);
void calculate_lambda(fluid_particle_t *fluid_particles, neighbors_t *neighbors_grid, param_t *params);
void update_dp(fluid_particle_t *fluid_particles, neighbors_t *neighbors_grid, param_t *params);
void identify_oob_particles(fluid_particle_t *fluid_particles, communication_t *communication, param_t *params);
void predict_positions(fluid_particle_t *fluid_particles, AABB_t *boundary_global, param_t *params);
void check_velocity(double *v_x, double *v_y, double *v_z);
void update_velocities(fluid_particle_t *fluid_particles, param_t *params);
void boundary_conditions(fluid_particle_t *fluid_particles, unsigned int i, AABB_t *boudnary);
void allocate_fluid(fluid_particle_t **fluid_particles, param_t *params);
void init_particles(fluid_particle_t *fluid_particles,
                neighbors_t *neighbors, AABB_t* water, AABB_t* boundary_global,
                param_t* params);


#endif
