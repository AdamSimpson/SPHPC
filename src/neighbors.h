#ifndef SPH_SRC_NEIGHBORS_H_
#define SPH_SRC_NEIGHBORS_H_

#include <stdbool.h>

// Forward Declaration
struct Particles;
struct Params;
struct AABB;

struct Neighbors {
  unsigned int *start_indices; // Start index for hash values
  unsigned int *end_indices;   // End index for hash values
  unsigned int *hash_values;   // Array of hash values
  unsigned int *particle_ids;  // Array of particle id's
  int max_neighbors;
  int hash_size_x;
  int hash_size_y;
  int hash_size_z;
  struct Neighbor *restrict particle_neighbors;
  double hash_spacing;
};

struct Neighbor {
    unsigned int neighbor_indices[60];
    int number_fluid_neighbors;
};

// These functions are to be visible to both C and C++
#ifdef __cplusplus
extern "C" {
#endif

void AllocateNeighbors(struct Neighbors *const neighbors,
                       const struct Particles *particles,
                       const struct Params *const params,
                       const struct AABB *const boundary_global);

void FreeNeighbors(struct Neighbors *neighbors);

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
              const struct Neighbors *const neighbors);

void FindCellBounds(const struct Particles *particles,
                    struct Neighbors *const neighbors);

void FillParticleNeighbors(struct Neighbors *const neighbors,
                           const struct Particles *particles,
                           const struct Params *const params,
                           const unsigned int p_index);

void FillNeighbors(const struct Particles *particles,
                   const struct Params *const params,
                   struct Neighbors *const neighbors);

#endif
