#include "neighbors.h"
#include "debug.h"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

void AllocateNeighbors(Neighbors *const neighbors,
                       const Params *const params,
                       const AABB *const boundary_global) {

  // Allocate neighbors array
  neighbors->particle_neighbors = calloc(params->max_fluid_particles_local,
                                         sizeof(Neighbor));
  if(neighbors->particle_neighbors == NULL)
    printf("Could not allocate neighbors\n");

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

  neighbors->hash = calloc(hash_size, sizeof(HashBucket));
  if(neighbors->hash == NULL)
    printf("Could not allocate hash\n");
  neighbors->start_indices = calloc(hash_size, sizeof(uint));
  if(neighbors->start_indices == NULL)
    printf("Could not allocate start_indices\n");
  neighbors->end_indices = calloc(hash_size, sizeof(uint));
  if(neighbors->end_indices == NULL)
    printf("Could not allocate end_indices\n");
  neighbors->hash_values = calloc(params->max_fluid_particles_local, sizeof(uint));
  if(neighbors->hash_values == NULL)
    printf("Could not allocate hash_values\n");
  neighbors->particle_ids = calloc(params->max_fluid_particles_local, sizeof(uint));
  if(neighbors->particle_ids == NULL)
    printf("Could not allocate particle_ids\n");
}

void FreeNeighbors(Neighbors *neighbors) {
  free(neighbors->particle_neighbors);
  free(neighbors->hash);
}

// Uniform grid hash
// We don't check if the position is out of bounds so x,y,z must be valid
uint HashVal(const Neighbors *const neighbors,
                     const double x,
                     const double y,
                     const double z) {
  const double spacing = neighbors->hash_spacing;

  // Calculate grid coordinates
  const uint grid_x = floor(x/spacing);
  const uint grid_y = floor(y/spacing);
  const uint grid_z = floor(z/spacing);

  // If using glboal boundary size this can be static
  const uint num_x = neighbors->hash_size_x;
  const uint num_y = neighbors->hash_size_y;

  const uint hash_val = (num_x * num_y * grid_z) + (grid_y * num_x + grid_x);

  return hash_val;
}

// Calculate and fill all neighbor particles
void FindAllNeighbors(fluid_sim_t *fluid_sim)
{
  HashParticles(fluid_sim);
  SortHash(fluid_sim);
  FindCellBounds(fluid_sim);
  FillNeighbors(fluid_sim);
}

// Hash all particles
void HashParticles(const Params *const params
                   const Neighbors *const neighbors) {

  uint *const hash_values = neighbors->hash_values;
  uint *const particle_ids = neighbors->particle_ids;
  int num_particles = params->number_fluid_particles_local + params->number_halo_particles;

  for (int i=0; i<num_particles; i++) {
    hash_values[i] =  HashVal(neighbors,
                              fluid_particles->x[p_index],
                              fluid_particles->y[p_index],
                              fluid_particles->z[p_index]);
    particle_ids[i] = i;
  }
}

// Sort list of particle id's based upon what their hash value is
void SortHash(const Params *const params
              const Neighbors *const neighbors) {

  uint *const keys = neighbrs->hash_values;
  uint *const values = neighbors->particle_ids;
  const int total_particles = params->number_fluid_particles_local + params->number_halo_particles;
  thrust::sort_by_key(thrust::host, keys, keys+total_particles, values);
}

// Find start and end of hash cells
// Method taken from NVIDIA SDK:
// http://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf
// Note that end index is one past the "end"
void FindCellBounds(const Params *const params
                    const Neighbors *const neighbors) {
  // Reset start indicies
  const int length_hash = neighbors->hash_size_x
                        * neighbors->hash_size_y
                        * neighbors->hash_size_z;
  memset(neighbors->start_indices, ((uint)-1), length_hash*sizeof(uint));

  const int num_particles = params->number_fluid_particles_local
                          + params->number_halo_particles;

  // If this particle has a different cell index to the previous
  // particle then it must be the first particle in the cell,
  // so store the index of this particle in the cell.
  // As it isn't the first particle, it must also be the cell end of
  // the previous particle's cell
  for (int i=0; i<num_particles; i++) {
    const uint hash = neighbors->hash_values[i];
    const uint i_prev = (i-1)>=0?(i-1):0;
    const uint hash_prev = neighbors->hash_values[i_prev];

    if (i==0 || hash!=hash_prev) {
      neighbors->start_indices[hash] = i;
      if (i > 0)
        neighbors->end_indices[hash_prev] = i;
    }
    if (i == num_particles - 1)
      neighbors->end_indices[hash] = i+1;
  }
}

// Neighbors are accessed multiple times per step so we keep them in buckets
void FillParticleNeighbors(const Neighbors *const neighbors,
                           const Params *const params,
                           const FluidParticles *particles,
                           const uint int p_index) {

  const int max_neighbors = neighbors->max_neighbors;
  const double smoothing_radius2 = params->smoothing_radius
                                 * params->smoothing_radius;
  const double spacing = neighbors->spacing;

  // Get neighbor bucket for particle p
  neighbor_t *const neighbors = neighbors[p_index];
  neighbors->number_fluid_neighbors = 0;

  // Calculate coordinates within bucket grid
  const uint grid_x = floor(particles->x[p_index]/spacing);
  const uint grid_y = floor(particles->y[p_index]/spacing);
  const uint grid_z = floor(particles->z[p_index]/spacing);

  // Go through neighboring grid buckets
  for (int dz=-1; dz<=1; ++dz) {
    for (int dy=-1; dy<=1; ++dy) {
      for (int dx=-1; dx<=1; ++dx) {
        // If the neighbor grid bucket is outside of the grid we don't process it
        if (grid_y+dy < 0 || grid_x+dx < 0 || grid_z+dz < 0
             || (grid_x+dx) >= neighbors->hash_size_x
             || (grid_y+dy) >= neighbors->hash_size_y
             || (grid_z+dz) >= neighbors->hash_size_z)
               continue;

         // Linear hash index for grid bucket
         const uint bucket_index = (grid_z+dz)
                                   *(neighbors->hash_size_x + neighbors->hash_size_y)
                                   +(grid_y+dy) *neighbors->hash_size_x + grid_x+dx;

         // Start index for hash value of current neighbor grid bucket
         uint start_index = neighbors->start_indices[bucket_index];

         // If neighbor grid bucket is not empty
         if (start_index != (uint)-1) {
           end_index = neighbors->end_indices[bucket_index];

           for (int j=start_index; j<end_index; ++j) {
             q_index = neighbors->particle_ids[j];

             // Continue if same particle
             if (p_index==q_index)
               continue;

             // Calculate distance squared
             const double r2 = (particles->x[p_index] - particles->x[q_index])
                              *(particles->x[p_index] - particles->x[q_index])
                              +(particles->y[p_index] - particles->y[q_index])
                              *(particles->y[p_index] - particles->y[q_index])
                              +(particles->z[p_index] - particles->z[q_index])
                              *(particles->z[p_index] - particles->z[q_index]);

             // If inside smoothing radius and enough space in p's neighbor bucket add q
             if(r2<smoothing_radius2 && neighbors->number_fluid_neighbors < max_neighbors)
                neighbors->neighbor_indices[neighbors->number_fluid_neighbors++] = q_index;
          }
        }
      } //dx
    } //dy
  } //dz
}

void fill_neighbors(const Neighbors *const neighbors,
                    const Params *const params,
                    const FluidParticles *particles) {
    const int total_particles = params->number_fluid_particles_local;

    // Fill neighbor bucket for all resident particles
    for (int i=0; i<total_particles; ++i) {
        FillParticleNeighbors(neighbors,
                              params,
                              particles,
                              i);
    }
}
