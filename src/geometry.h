#ifndef SPH_SRC_GEOMETRY_H_
#define SPH_SRC_GEOMETRY_H_

typedef struct AABB AABB_t;

#include "fluid.h"
#include "communication.h"
#include "simulation.h"

struct AABB {
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double min_z;
    double max_z;
};//Axis aligned bounding box

double min(double a, double b);
double max(double a, double b);
int sgn(double x);
void construct_fluid_volume(fluid_particle_t *fluid_particles, AABB_t* fluid,
                            param_t *params);
void set_particle_numbers(AABB_t *fluid_global,
                          communication_t *communications, param_t *params);
void partition_problem(AABB_t *boundary_global, AABB_t *fluid_global,
                       int *x_start, int *number_particles_x, param_t *params);
void check_partition(fluid_particle_t *fluid_particles,
                     param_t *params);
#endif
