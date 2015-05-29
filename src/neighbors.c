#include "neighbors.h"
#include "debug.h"
#include "safe_alloc.h"
#include "particles.h"
#include "geometry.h"
#include "simulation.h"
#include "thrust_c.h"
#include <math.h>
#include <stdbool.h>
#include <string.h>

void AllocInitNeighbors(struct Neighbors *const neighbors,
                        const struct Particles *particles,
                        const struct Params *const params,
                        const struct AABB *const boundary_global) {

  // Allocate neighbors array
  neighbors->neighbor_buckets = (struct NeighborBucket*)
                                  SAFE_ALLOC(particles->max_local,
                                             sizeof(struct NeighborBucket));

  // +1 added because range begins at 0
  neighbors->hash_size_x = ceil((boundary_global->max_x - boundary_global->min_x)
                         / params->smoothing_radius) + 1;
  neighbors->hash_size_y = ceil((boundary_global->max_y - boundary_global->min_y)
                         / params->smoothing_radius) + 1;
  neighbors->hash_size_z = ceil((boundary_global->max_z - boundary_global->min_z)
                         / params->smoothing_radius) + 1;
  size_t hash_size = neighbors->hash_size_x
                   * neighbors->hash_size_y
                   * neighbors->hash_size_z;

  neighbors->start_indices = (unsigned int*)SAFE_ALLOC(hash_size,
                                                       sizeof(unsigned int));
  neighbors->end_indices = (unsigned int*)SAFE_ALLOC(hash_size,
                                                     sizeof(unsigned int));
  neighbors->hash_values = (unsigned int*)SAFE_ALLOC(particles->max_local,
                                                     sizeof(unsigned int));
  neighbors->particle_ids = (unsigned int*)SAFE_ALLOC(particles->max_local,
                                                      sizeof(unsigned int));
}

void FinalizeNeighbors(struct Neighbors *neighbors) {
  free(neighbors->neighbor_buckets);
  free(neighbors->start_indices);
  free(neighbors->end_indices);
  free(neighbors->hash_values);
  free(neighbors->particle_ids);
}

// Calculate and fill all neighbor particles
void FindAllNeighbors(const struct Particles *const particles,
                      const struct Params *const params,
                      struct Neighbors *const neighbors) {
  HashParticles(particles, neighbors);
  SortHash(particles, params, neighbors);
  FindCellBounds(particles, neighbors);
  FillNeighbors(particles, params, neighbors);
}

// Uniform grid hash
// We don't check if the position is out of bounds so x,y,z must be valid
unsigned int HashVal(const struct Neighbors *const neighbors,
                     const double x,
                     const double y,
                     const double z) {
  const double spacing = neighbors->hash_spacing;

  // Calculate grid coordinates
  const unsigned int grid_x = floor(x/spacing);
  const unsigned int grid_y = floor(y/spacing);
  const unsigned int grid_z = floor(z/spacing);

  // If using glboal boundary size this can be static
  const unsigned int num_x = neighbors->hash_size_x;
  const unsigned int num_y = neighbors->hash_size_y;

  const unsigned int hash_val = (num_x * num_y * grid_z)
                              + (grid_y * num_x + grid_x);

  return hash_val;
}

// Hash all particles
void HashParticles(const struct Particles *const particles,
                   struct Neighbors *const neighbors) {

  unsigned int *const hash_values = neighbors->hash_values;
  unsigned int *const particle_ids = neighbors->particle_ids;
  const int num_particles = particles->local_count
                          + particles->halo_count_left
                          + particles->halo_count_right;

  for (int i=0; i<num_particles; i++) {
    hash_values[i] =  HashVal(neighbors,
                              particles->x_star[i],
                              particles->y_star[i],
                              particles->z_star[i]);
    particle_ids[i] = i;
  }
}

// Sort list of particle id's based upon what their hash value is
void SortHash(const struct Particles *particles,
              const struct Params *const params,
              struct Neighbors *const neighbors) {

  unsigned int *const keys = neighbors->hash_values;
  unsigned int *const values = neighbors->particle_ids;
  const int total_particles = particles->local_count
                            + particles->halo_count_left
                            + particles->halo_count_right;
  SortByKey(keys, values, total_particles);

}

// Find start and end of hash cells
// Method taken from NVIDIA SDK:
// http://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf
// Note that end index is one past the "end"
void FindCellBounds(const struct Particles *particles,
                    struct Neighbors *const neighbors) {
  // Reset start indicies
  const int length_hash = neighbors->hash_size_x
                        * neighbors->hash_size_y
                        * neighbors->hash_size_z;

  const int num_particles = particles->local_count
                          + particles->halo_count_left
                          + particles->halo_count_right;

  // Find start and end indices for each
  FindLowerBounds(neighbors->hash_values,
                  num_particles,
                  length_hash,
                  neighbors->start_indices);
  FindUpperBounds(neighbors->hash_values,
                  num_particles,
                  length_hash,
                  neighbors->end_indices);
}

// Neighbors are accessed multiple times per step so we keep them in buckets
void FillParticleNeighbors(struct Neighbors *const neighbors,
                           const struct Params *const params,
                           const struct Particles *particles,
                           const unsigned int p_index) {

  const int max_neighbors = neighbors->max_neighbors;
  const double smoothing_radius2 = params->smoothing_radius
                                 * params->smoothing_radius;
  const double spacing = neighbors->hash_spacing;

  // Get neighbor bucket for particle p
  struct NeighborBucket *const neighbor_bucket = &neighbors->neighbor_buckets[p_index];
  neighbor_bucket->count = 0;

  const double px = particles->x_star[p_index];
  const double py = particles->y_star[p_index];
  const double pz = particles->z_star[p_index];

  // Go through neighboring grid buckets
  for (int dz=-1; dz<=1; ++dz) {
    const double z = pz + dz*spacing;
    for (int dy=-1; dy<=1; ++dy) {
      const double y = py + dy*spacing;
      for (int dx=-1; dx<=1; ++dx) {
        const double x = px + dx*spacing;

        // Make sure that the position is valid
        if (floor(x/spacing) > neighbors->hash_size_x-1 || x < 0.0 ||
            floor(y/spacing) > neighbors->hash_size_y-1 || y < 0.0 ||
            floor(z/spacing) > neighbors->hash_size_z-1 || z < 0.0)
          continue;

        // Calculate hash index at neighbor point
        const int neighbor_index = HashVal(neighbors, x, y, z);

        // Start index for hash value of current neighbor grid bucket
        const unsigned int start_index = neighbors->start_indices[neighbor_index];
        const unsigned int end_index = neighbors->end_indices[neighbor_index];
        // If neighbor grid bucket is not empty
        if (start_index != end_index) {

          for (int j=start_index; j<end_index; ++j) {
            const int q_index = neighbors->particle_ids[j];

            // Continue if same particle
            if (p_index==q_index)
              continue;

            // Calculate distance squared
            const double x_diff = particles->x_star[p_index]
                                - particles->x_star[q_index];
            const double y_diff = particles->y_star[p_index]
                                - particles->y_star[q_index];
            const double z_diff = particles->z_star[p_index]
                                - particles->z_star[q_index];
            const double r2 = x_diff*x_diff + y_diff*y_diff + z_diff*z_diff;

            // If inside smoothing radius and enough space
            // in p's neighbor bucket add q
            if(r2<smoothing_radius2 &&
               neighbor_bucket->count < max_neighbors) {
                const int num_neighbors = neighbor_bucket->count;
                neighbor_bucket->neighbor_indices[num_neighbors] = q_index;
                ++(neighbor_bucket->count);
            }
          }
        }
      } //dx
    } //dy
  } //dz
}

void FillNeighbors(const struct Particles *particles,
                   const struct Params *const params,
                   struct Neighbors *const neighbors) {
  const int total_particles = particles->local_count;

  // Fill neighbor bucket for all resident particles
  for (int i=0; i<total_particles; ++i) {
    FillParticleNeighbors(neighbors,
                          params,
                          particles,
                          i);
  }
}
