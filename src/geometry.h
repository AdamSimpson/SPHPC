#ifndef SPH_SRC_GEOMETRY_H_
#define SPH_SRC_GEOMETRY_H_

typedef struct AABB AABB;

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

void construct_fluid_volume(FluidParticle *const fluid_particles,
                            Params *const params,
                            const AABB *const fluid);

void set_particle_numbers(const AABB *const fluid_global,
                          Params *const params,
                          Communication *const communication);

void check_partition(Params *const params);

#endif
