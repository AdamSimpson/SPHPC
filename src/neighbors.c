#include "neighbors.h"

#include <math.h>
#include <stdbool.h>
#include <string.h>

void AllocateNeighbors(Neighbors *const neighbors,
                        const Params *const params,
                        const AABB *const boundary_global)
{
    // Allocate neighbors array
    neighbors->particle_neighbors = calloc(params->max_fluid_particles_local, sizeof(Neighbor));
    if(neighbors->particle_neighbors == NULL)
        printf("Could not allocate neighbors\n");

    // +1 added because range begins at 0
    neighbors->hash_size_x = ceil((boundary_global->max_x - boundary_global->min_x) / params->smoothing_radius) + 1;
    neighbors->hash_size_y = ceil((boundary_global->max_y - boundary_global->min_y) / params->smoothing_radius) + 1;
    neighbors->hash_size_z = ceil((boundary_global->max_z - boundary_global->min_z) / params->smoothing_radius) + 1;
    size_t hash_size = neighbors->hash_size_x * neighbors->hash_size_y * neighbors->hash_size_z;

    neighbors->hash = calloc(hash_size, sizeof(HashBucket));
    if(neighbors->hash == NULL)
        printf("Could not allocate hash\n");
}

void FreeNeighbors(Neighbors *neighbors)
{
  free(neighbors->particle_neighbors);
  free(neighbors->hash);
}

// Uniform grid hash
// We don't check if the position is out of bounds so x,y,z must be valid
unsigned int HashVal(const Neighbors *const neighbors,
                      const double x,
                      const double y,
                      const double z)
{
    const double spacing = neighbors->hash_spacing;

    // Calculate grid coordinates
    const int grid_x = floor(x/spacing);
    const int grid_y = floor(y/spacing);
    const int grid_z = floor(z/spacing);

    // If using glboal boundary size this can be static
    const int num_x = neighbors->hash_size_x;
    const int num_y = neighbors->hash_size_y;

    const int grid_position = (num_x * num_y * grid_z) + (grid_y * num_x + grid_x);

    return grid_position;
}

// Add halo particles to neighbors array
// Halo particles are not added to hash, only neighbors list
// Neighbors may be more than h away...since distance is computed in all smoothing functions
// it is a waste to check as we hash as well
void HashHalo(const FluidParticle *const fluid_particles,
               const Params *const params,
               const AABB *const boundary,
               Neighbors *const neighbors)
{
    const int n_s = params->number_fluid_particles_local;
    const int n_f = n_s + params->number_halo_particles_left + params->number_halo_particles_right;
    const double spacing = neighbors->hash_spacing;
    const double h = params->smoothing_radius;

    for(int i=n_s; i<n_f; i++)
    {
        const FluidParticle *const h_p = &fluid_particles[i];

        // Search bins around current particle
        for (int dx=-1; dx<=1; dx++) {
            const double x = h_p->x_star + dx*spacing;
            for (int dy=-1; dy<=1; dy++) {
                const double y = h_p->y_star + dy*spacing;
                for (int dz=-1; dz<=1; dz++) {
                    const double z = h_p->z_star + dz*spacing;

                    // Make sure that the position is valid
                    if( floor(x/spacing) > neighbors->hash_size_x-1 || x < 0 ||
                        floor(y/spacing) > neighbors->hash_size_y-1 || y < 0 ||
                        floor(z/spacing) > neighbors->hash_size_z-1 || z < 0)
                      continue;

                    // Calculate hash index at neighbor point
                    const int index = HashVal(neighbors, x,y,z);
                      // Go through each fluid particle in neighbor point bucket
                      for (int n=0;n<neighbors->hash[index].number_fluid;n++) {
                          const FluidParticle *const q = neighbors->hash[index].fluid_particles[n];
                          const double r = sqrt((h_p->x_star-q->x_star)*(h_p->x_star-q->x_star)
                                 + (h_p->y_star-q->y_star)*(h_p->y_star-q->y_star)
                                 + (h_p->z_star-q->z_star)*(h_p->z_star-q->z_star));
                          if(r > h)
                              continue;

                          // Get neighbor ne for particle q
                          Neighbor *const ne = &neighbors->particle_neighbors[q->id];
                          // Make sure not to add duplicate neighbors
                          bool duped = false;
                          for (int dupes=0; dupes < ne->number_fluid_neighbors; dupes++) {
                                if (ne->neighbor_indices[dupes] == i) {
                                    duped = true;
                                    break;
                                }
                            }
                            if (!duped && ne->number_fluid_neighbors < 200) {
                                ne->neighbor_indices[ne->number_fluid_neighbors] = i;
                                ne->number_fluid_neighbors++;
                            }

                      }
                }
            }
        }
    }
}

// Fill fluid particles into hash
// Neighbors may be more than h away...since distance is computed in all smoothing functions
// it is a waste to check as we hash as well
void HashFluid(const FluidParticle *const fluid_particles,
                const Params *const params,
                const AABB *const boundary,
                Neighbors *neighbors)
{
//        int i,dx,dy,dz,n,c;
//        double x,y,z, px,py,pz;
        const double spacing = neighbors->hash_spacing;
        const double h = params->smoothing_radius;
        const int n_f = params->number_fluid_particles_local;
//        const FluidParticle *p, *q, *q_neighbor;
//        Neighbor *ne;
//        double r;
//        unsigned int index, neighbor_index;

        const unsigned int length_hash = neighbors->hash_size_x
                                 * neighbors->hash_size_y
                                 * neighbors->hash_size_z;

        // zero out number of particles in bucket
        for (int index=0; index<length_hash; index++){
            neighbors->hash[index].number_fluid = 0;
            neighbors->hash[index].hashed = false;
        }

        // First pass - insert fluid particles into hash
        for (int i=0; i<n_f; i++) {
            const FluidParticle *const p = &fluid_particles[i];

            neighbors->particle_neighbors[i].number_fluid_neighbors = 0;

            const int index = HashVal(neighbors, p->x_star, p->y_star, p->z_star);

            if (neighbors->hash[index].number_fluid < neighbors->max_neighbors) {
                neighbors->hash[index].fluid_particles[neighbors->hash[index].number_fluid] = p;
                neighbors->hash[index].number_fluid++;
            }
        }

        // Second pass - fill particle neighbors by processing buckets
      	// Could also iterate through hash directly but particles array will be shorter
        for (int i=0; i<n_f; i++) {
           // Calculate hash index of bucket
           const FluidParticle *const p = &fluid_particles[i];
           const double px = p->x_star;
           const double py = p->y_star;
           const double pz = p->z_star;

           const int index = HashVal(neighbors, px, py, pz);

           // If this bucket has been done try the next one
           if(neighbors->hash[index].hashed)
               continue;

            // Check neighbors of current bucket
            for (int dx=-1; dx<=1; dx++) {
                const double x = px + dx*spacing;
                for (int dy=-1; dy<=1; dy++) {
                    const double y = py + dy*spacing;
                    for (int dz=-1; dz<=1; dz++) {
                        const double z = pz + dz*spacing;

                        // Make sure that the position is valid
                        if( floor(x/spacing) > neighbors->hash_size_x-1 || x < 0 ||
                            floor(y/spacing) > neighbors->hash_size_y-1 || y < 0 ||
                            floor(z/spacing) > neighbors->hash_size_z-1 || z < 0)
                          continue;

                        // Calculate hash index at neighbor point
                        const int neighbor_index = HashVal(neighbors, x, y, z);

                        // Add neighbor particles to each particle in current bucket
                        for (int c=0; c<neighbors->hash[index].number_fluid; c++) {
			                      // Particle in currently being worked on bucket
                            const FluidParticle *const q = neighbors->hash[index].fluid_particles[c];
                            Neighbor *const ne = &neighbors->particle_neighbors[q->id];
			                      for(int n=0; n<neighbors->hash[neighbor_index].number_fluid; n++){
                                // Append neighbor to q's neighbor list
		   	                        const FluidParticle *const q_neighbor = neighbors->hash[neighbor_index].fluid_particles[n];
                                if(q->id == q_neighbor->id)
                                     continue;

                                const double r = sqrt((q_neighbor->x_star-q->x_star)*(q_neighbor->x_star-q->x_star)
                                       + (q_neighbor->y_star-q->y_star)*(q_neighbor->y_star-q->y_star)
                                       + (q_neighbor->z_star-q->z_star)*(q_neighbor->z_star-q->z_star));
                                if(r > h)
                                    continue;

                                if(ne->number_fluid_neighbors > neighbors->max_neighbors)
                                    printf("too many neighbors: %f, %f, %f\n", fluid_particles[i].x_star,fluid_particles[i].y_star,fluid_particles[i].z_star);
                                    ne->neighbor_indices[ne->number_fluid_neighbors++] = q_neighbor->id;
		                        }
                       }

                   } // end dz
                } // end dy
             }  // end dx

	    // This bucket has been hashed and does not need hashed again
	    neighbors->hash[index].hashed = true;

      } // end main particle loop

}// end function
