#ifndef NEIGHBORS_h
#define NEIGHBORS_h

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
    fluid_particle_t *fluid_particles[200];
    unsigned int number_fluid;
    bool hashed;
};

unsigned int hash_val(neighbors_t *neighbors, double x, double y, double z);
void hash_fluid(fluid_particle_t *fluid_particles, neighbors_t *neighbors, AABB_t *boundary, param_t *params);
void hash_halo(fluid_particle_t *fluid_particles, neighbors_t *neighbors, AABB_t *boundary, param_t *params);
void allocate_hash(neighbors_t *neighbors, AABB_t *boundary_global, param_t *params);

#endif
