#include "neighbors.h"
#include "debug.h"
#include "safe_alloc.h"
#include "particles.h"
#include "geometry.h"
#include "simulation.h"
#include "thrust_c.h"
#include "openacc.h"
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

  neighbors->hash_spacing = params->smoothing_radius;

  #pragma acc enter data copyin( neighbors[:1],                               \
                         neighbors->neighbor_buckets[0:particles->max_local], \
                         neighbors->start_indices[0:hash_size],               \
                         neighbors->end_indices[0:hash_size],                 \
                         neighbors->hash_values[0:particles->max_local],      \
                         neighbors->particle_ids[0:particles->max_local])
}

void FinalizeNeighbors(struct Neighbors *neighbors) {

  #pragma acc exit data delete( neighbors->neighbor_buckets, \
    neighbors->start_indices,                                \
    neighbors->end_indices,                                  \
    neighbors->hash_values,                                  \
    neighbors->particle_ids)

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
#pragma acc routine seq
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

  // OpenACC can't reliably handle SoA...
  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;

  #pragma acc parallel loop present(x_star, y_star, z_star, \
                                    neighbors, hash_values, particle_ids) default(none)
  for (int i=0; i<num_particles; ++i) {
    hash_values[i] =  HashVal(neighbors,
                              x_star[i],
                              y_star[i],
                              z_star[i]);
    particle_ids[i] = i;
  }
}

// Sort list of particle id's based upon what their hash value is
void SortHash(const struct Particles *particles,
              const struct Params *const params,
              struct Neighbors *const neighbors) {

  unsigned int *const hash_values = neighbors->hash_values;
  unsigned int *const particle_ids = neighbors->particle_ids;
  const int local_count = particles->local_count
                            + particles->halo_count_left
                            + particles->halo_count_right;

  // OpenACC can't reliably handle SoA...
  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;
  double *x    = particles->x;
  double *y    = particles->y;
  double *z    = particles->z;
  double *v_x    = particles->v_x;
  double *v_y    = particles->v_y;
  double *v_z    = particles->v_z;
  int *id         = particles->id;

  #pragma acc host_data use_device(hash_values, particle_ids)
  {
    void *cuda_stream = acc_get_cuda_stream(acc_async_sync);

    SortByKey(hash_values, particle_ids, local_count, cuda_stream);
  }
/*
  // This is broke with -O3...
  // reorder particles in particle_id order for better memory access
  // This will fuckup halo indices
  #pragma acc host_data use_device(hash_values, x_star, y_star, z_star,  \
                                   x, y, z, v_x, v_y, v_z)
  {
    void *cuda_stream = acc_get_cuda_stream(acc_async_sync);

  SortParticlesByKey(hash_values,
                     local_count,
                     x_star, y_star, z_star,
                     x, y, z,
                     v_x, v_y, v_z, cuda_stream);
  }

  // Reset particle id's
  #pragma acc parallel loop present(id)
  for(int i=0; i<local_count; i++) {
    id[i] = i;
  }
*/
}

// Find start and end of hash cells
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

  // OpenACC doesn't accept struct member pointers in use_device
  unsigned int *hash_values = neighbors->hash_values;
  unsigned int *start_indices = neighbors->start_indices;
  unsigned int *end_indices = neighbors->end_indices;

  // Find start and end indices for each
  #pragma acc host_data use_device(hash_values, start_indices, end_indices)
  {
        void *cuda_stream = acc_get_cuda_stream(acc_async_sync);

    FindLowerBounds(hash_values,
                    num_particles,
                    length_hash,
                    start_indices,
                    cuda_stream);

    FindUpperBounds(hash_values,
                    num_particles,
                    length_hash,
                    end_indices,
                    cuda_stream);
  }
}

void FillNeighbors(const struct Particles *particles,
                   const struct Params *const params,
                   struct Neighbors *const neighbors) {

  const int num_particles = particles->local_count;

  const double smoothing_radius2 = params->smoothing_radius
                                 * params->smoothing_radius;
  const double spacing = neighbors->hash_spacing;

  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;
  struct NeighborBucket *neighbor_buckets = neighbors->neighbor_buckets;

  // Fill neighbor bucket for all resident particles
  #pragma acc parallel loop present(x_star, y_star, z_star, neighbors, neighbor_buckets)
  for (int i=0; i<num_particles; ++i) {
      const int p_index = i;

      // Get neighbor bucket for particle p
      struct NeighborBucket *const neighbor_bucket = &neighbor_buckets[p_index];
      neighbor_bucket->count = 0;

      const double px = x_star[p_index];
      const double py = y_star[p_index];
      const double pz = z_star[p_index];

      // Go through neighboring grid buckets
      #pragma acc loop seq collapse(3)
      for (int dz=-1; dz<=1; ++dz) {
        for (int dy=-1; dy<=1; ++dy) {
          for (int dx=-1; dx<=1; ++dx) {

            const double y = py + dy*spacing;
            const double z = pz + dz*spacing;
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
                const double x_diff = px
                                    -x_star[q_index];
                const double y_diff = py
                                    -y_star[q_index];
                const double z_diff = pz
                                    -z_star[q_index];
                const double r2 = x_diff*x_diff + y_diff*y_diff + z_diff*z_diff;

                // If inside smoothing radius and enough space
                // in p's neighbor bucket add q
                if(r2<smoothing_radius2 &&
                   neighbor_bucket->count < MAX_NEIGHBORS) {
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

}
