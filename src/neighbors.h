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

void FindAllNeighbors(const struct Particles *const particles,
                      const struct Params *const params,
                      struct Neighbors *const neighbors);
#ifdef __cplusplus
}
#endif

unsigned int HashVal(const struct Neighbors *const neighbors,
                     const double x,
                     const double y,
                     const double z);

void HashParticles(const struct Particles *const particles,
                   struct Neighbors *const neighbors);

void SortHash(const struct Particles *particles,
              const struct Params *const params,
              struct Neighbors *const neighbors);

void FindCellBounds(const struct Particles *particles,
                    struct Neighbors *const neighbors);

void FillParticleNeighbors(struct Neighbors *const neighbors,
                           const struct Params *const params,
                           const struct Particles *particles,
                           const unsigned int p_index);

void FillNeighbors(const struct Particles *particles,
                   const struct Params *const params,
                   struct Neighbors *const neighbors);

#endif
