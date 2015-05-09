#ifndef SPH_SRC_NEIGHBORS_H_
#define SPH_SRC_NEIGHBORS_H_

typedef struct BUCKET bucket_t;
typedef struct NEIGHBORS neighbors_t;
typedef struct NEIGHBOR neighbor_t;

#include <stdbool.h>

#include "fluid.h"
#include "geometry.h"
#include "simulation.h"

struct NEIGHBORS {
  neighbor_t *particle_neighbors;
  bucket_t *hash;
  double hash_spacing;
  int max_neighbors;
  int hash_size_x;
  int hash_size_y;
  int hash_size_z;
};

struct NEIGHBOR {
    unsigned int neighbor_indices[60];
    int number_fluid_neighbors;
};

// hash 'bucket' for hash value
struct BUCKET {
    const fluid_particle_t *fluid_particles[200];
    unsigned int number_fluid;
    bool hashed;
};

void allocate_neighbors(neighbors_t *const neighbors,
                        const param_t *const params,
                        const AABB_t *const boundary_global);

void free_neighbors(neighbors_t *neighbors);

unsigned int hash_val(const neighbors_t *const neighbors,
                      const double x,
                      const double y,
                      const double z);

void hash_halo(const fluid_particle_t *const fluid_particles,
               const param_t *const params,
               const AABB_t *const boundary,
               neighbors_t *const neighbors);

void hash_fluid(const fluid_particle_t *const fluid_particles,
                const param_t *const params,
                const AABB_t *const boundary,
                neighbors_t *neighbors);

#endif
