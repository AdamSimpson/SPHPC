#ifndef SPH_SRC_NEIGHBORS_H_
#define SPH_SRC_NEIGHBORS_H_

#include <stdbool.h>

#define MAX_NEIGHBORS 60

// Forward Declarations
struct Particles;
struct Params;
struct AABB;

struct Neighbors {
  unsigned int *restrict start_indices; // Start index for hash values
  unsigned int *restrict end_indices;   // End index for hash values
  unsigned int *restrict hash_values;   // Array of hash values
  unsigned int *restrict particle_ids;  // Array of particle id's
  int hash_size_x;
  int hash_size_y;
  int hash_size_z;
  struct NeighborBucket *restrict neighbor_buckets;
  double hash_spacing;
};

struct NeighborBucket {
    unsigned int neighbor_indices[MAX_NEIGHBORS];
    int count;
};

// These functions are to be visible to both C and C++
#ifdef __cplusplus
extern "C" {
#endif

void AllocInitNeighbors(struct Neighbors *const neighbors,
                       const struct Particles *particles,
                       const struct Params *const params,
                       const struct AABB *const boundary_global);

void FinalizeNeighbors(struct Neighbors *neighbors);

void FindAllNeighbors(const struct Particles *restrict particles,
                      const struct Params *restrict params,
                      struct Neighbors *restrict neighbors);
#ifdef __cplusplus
}
#endif

unsigned int HashVal(const struct Neighbors *restrict neighbors,
                     const double x,
                     const double y,
                     const double z);

void HashParticles(const struct Particles *restrict particles,
                   struct Neighbors *restrict neighbors);

void SortHash(const struct Particles *restrict particles,
              const struct Params *restrict params,
              struct Neighbors *restrict neighbors);

void FindCellBounds(const struct Particles *restrict particles,
                    struct Neighbors *restrict neighbors);

void FillNeighbors(const struct Particles *restrict particles,
                   const struct Params *restrict params,
                   struct Neighbors *restrict neighbors);

#endif
