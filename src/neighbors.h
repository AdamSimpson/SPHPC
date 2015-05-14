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
  unsigned int *start_indices; // Start index for hash values
  unsigned int *end_indices;   // End index for hash values
  unsigned int *hash_values;   // Array of hash values
  unsigned int *particle_ids;  // Array of particle id's
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

// These functions are to be visible to both
#ifdef __cplusplus
extern "C" {
#endif
void AllocateNeighbors(Neighbors *const neighbors,
                       const Params *const params,
                       const AABB *const boundary_global);

void FreeNeighbors(Neighbors *neighbors);

void FindAllNeighbors(const Params *const params,
                      const FluidParticles *const particles,
                      Neighbors *const neighbors);
#ifdef __cplusplus
}
#endif

unsigned int HashVal(const Neighbors *const neighbors,
                     const double x,
                     const double y,
                     const double z);

void HashParticles(const Params *const params,
                   const FluidParticles *const particles,
                   Neighbors *const neighbors);

void SortHash(const Params *const params,
              const Neighbors *const neighbors);

void FindCellBounds(const Params *const params,
                    Neighbors *const neighbors);

void FillParticleNeighbors(Neighbors *const neighbors,
                           const Params *const params,
                           const FluidParticles *particles,
                           const unsigned int p_index);

void FillNeighbors(const Params *const params,
                   const FluidParticles *particles,
                   Neighbors *const neighbors);

#endif
