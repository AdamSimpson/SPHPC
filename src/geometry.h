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

void construct_fluid_volume(fluid_particle_t *const fluid_particles,
                            param_t *const params,
                            const AABB_t *const fluid);

void set_particle_numbers(const AABB_t *const fluid_global,
                          communication_t *const communication,
                          param_t *const params);

void check_partition(param_t *const params);

#endif
