#ifndef SPH_SRC_NEIGHBORS_H_
#define SPH_SRC_NEIGHBORS_H_

typedef struct BUCKET HashBucket;
typedef struct NEIGHBORS Neighbors;
typedef struct NEIGHBOR Neighbor;

#include <stdbool.h>

#include "fluid.h"
#include "geometry.h"
#include "simulation.h"

struct NEIGHBORS {
  Neighbor *particle_neighbors;
  HashBucket *hash;
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
    const FluidParticle *fluid_particles[200];
    unsigned int number_fluid;
    bool hashed;
};

void AllocateNeighbors(Neighbors *const neighbors,
                       const Params *const params,
                       const AABB *const boundary_global);

void FreeNeighbors(Neighbors *neighbors);

unsigned int HashVal(const Neighbors *const neighbors,
                     const double x,
                     const double y,
                     const double z);

void HashHalo(const FluidParticle *const fluid_particles,
              const Params *const params,
              const AABB *const boundary,
              Neighbors *const neighbors);

void HashFluid(const FluidParticle *const fluid_particles,
               const Params *const params,
               const AABB *const boundary,
               Neighbors *neighbors);

#endif
