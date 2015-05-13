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
  uint *start_indices; // Start index for hash values
  uint *end_indices;   // End index for hash values
  uint *hash_values;   // Array of hash values
  uint *particle_ids;  // Array of particle id's
  int max_neighbors;
  int hash_size_x;
  int hash_size_y;
  int hash_size_z;
  Neighbor *particle_neighbors;
  double hash_spacing;
};

struct NEIGHBOR {
    unsigned int neighbor_indices[60];
    int number_fluid_neighbors;
};

void AllocateNeighbors(Neighbors *const neighbors,
                       const Params *const params,
                       const AABB *const boundary_global);

void FreeNeighbors(Neighbors *neighbors);

unsigned int HashVal(const Neighbors *const neighbors,
                     const double x,
                     const double y,
                     const double z);

void HashHalo(const FluidParticles *const fluid_particles,
              const Params *const params,
              const AABB *const boundary,
              Neighbors *const neighbors);

void HashFluid(const FluidParticles *const fluid_particles,
               const Params *const params,
               const AABB *const boundary,
               Neighbors *neighbors);

#endif
