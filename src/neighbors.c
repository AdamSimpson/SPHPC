#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "neighbors.h"

void allocate_neighbors(neighbors_t *neighbors, AABB_t *boundary_global, param_t *params)
{
    // Allocate neighbors array
    neighbors->particle_neighbors = calloc(params->max_fluid_particles_local, sizeof(neighbor_t));
    if(neighbors->particle_neighbors == NULL)
        printf("Could not allocate neighbors\n");

    // +1 added because range begins at 0
    neighbors->hash_size_x = ceil((boundary_global->max_x - boundary_global->min_x) / params->smoothing_radius) + 1;
    neighbors->hash_size_y = ceil((boundary_global->max_y - boundary_global->min_y) / params->smoothing_radius) + 1;
    neighbors->hash_size_z = ceil((boundary_global->max_z - boundary_global->min_z) / params->smoothing_radius) + 1;
    unsigned int hash_size = neighbors->hash_size_x * neighbors->hash_size_y * neighbors->hash_size_z;

    neighbors->hash = calloc(hash_size, sizeof(bucket_t));
    if(neighbors->hash == NULL)
        printf("Could not allocate hash\n");
}

// Uniform grid hash
// We don't check if the position is out of bounds so x,y,z must be valid
inline unsigned int hash_val(neighbors_t *neighbors, double x, double y, double z)
{
    double spacing = neighbors->hash_spacing;

    // Calculate grid coordinates
    unsigned int grid_x,grid_y,grid_z;
    grid_x = floor(x/spacing);
    grid_y = floor(y/spacing);
    grid_z = floor(z/spacing);

    // If using glboal boundary size this can be static
    int num_x = neighbors->hash_size_x;
    int num_y = neighbors->hash_size_y;

    unsigned int grid_position = (num_x * num_y * grid_z) + (grid_y * num_x + grid_x);

    return grid_position;
}

// Add halo particles to neighbors array
// Halo particles are not added to hash, only neighbors list
// Neighbors may be more than h away...since distance is computed in all smoothing functions
// it is a waste to check as we hash as well
void hash_halo(fluid_particle_t *fluid_particles, neighbors_t *neighbors, AABB_t *boundary, param_t *params)
{
    printf("Enter hash halo\n");
    int index,i,dx,dy,dz,dupes,n;
    double x,y,z,r;
    bool duped;
    fluid_particle_t *h_p, *q;
    int n_s = params->number_fluid_particles_local;
    int n_f = n_s + params->number_halo_particles_left + params->number_halo_particles_right;
    double spacing = neighbors->hash_spacing;
    double h = params->smoothing_radius;
    neighbor_t *ne;

    for(i=n_s; i<n_f; i++)
    {
        h_p = &fluid_particles[i];

        // Search bins around current particle
        for (dx=-1; dx<=1; dx++) {
            x = h_p->x_star + dx*spacing;
            for (dy=-1; dy<=1; dy++) {
                y = h_p->y_star + dy*spacing;
                for (dz=-1; dz<=1; dz++) {
                    z = h_p->z_star + dz*spacing;

                    // Make sure that the position is valid
                    if( floor(x/spacing) > neighbors->hash_size_x-1 || x < 0 ||
                        floor(y/spacing) > neighbors->hash_size_y-1 || y < 0 ||
                        floor(z/spacing) > neighbors->hash_size_z-1 || z < 0)
                      continue;

                    // Calculate hash index at neighbor point
                    index = hash_val(neighbors, x,y,z);
                      // Go through each fluid particle in neighbor point bucket
                      for (n=0;n<neighbors->hash[index].number_fluid;n++) {
                          q = neighbors->hash[index].fluid_particles[n];
                          r = sqrt((h_p->x_star-q->x_star)*(h_p->x_star-q->x_star)
                                 + (h_p->y_star-q->y_star)*(h_p->y_star-q->y_star)
                                 + (h_p->z_star-q->z_star)*(h_p->z_star-q->z_star));
                          if(r > h)
                              continue;

                          // Get neighbor ne for particle q
                          ne = &neighbors->particle_neighbors[q->id];
                          // Make sure not to add duplicate neighbors
                          duped = false;
                          for (dupes=0; dupes < ne->number_fluid_neighbors; dupes++) {
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
void hash_fluid(fluid_particle_t *fluid_particles, neighbors_t *neighbors,
                AABB_t *boundary, param_t *params)
{
        printf("Hash fluid\n");
        int i,dx,dy,dz,n,c;
        double x,y,z, px,py,pz;
        double spacing = neighbors->hash_spacing;
        double h = params->smoothing_radius;
        int n_f = params->number_fluid_particles_local;
        fluid_particle_t *p, *q, *q_neighbor;
        neighbor_t *ne;
        double r;
        unsigned int index, neighbor_index;

        unsigned int length_hash = neighbors->hash_size_x
                                 * neighbors->hash_size_y
                                 * neighbors->hash_size_z;

        // zero out number of particles in bucket
        for (index=0; index<length_hash; index++){
            neighbors->hash[index].number_fluid = 0;
            neighbors->hash[index].hashed = false;
        }

        // First pass - insert fluid particles into hash
        for (i=0; i<n_f; i++) {
            p = &fluid_particles[i];

            neighbors->particle_neighbors[i].number_fluid_neighbors = 0;

            index = hash_val(neighbors, p->x_star, p->y_star, p->z_star);

            if (neighbors->hash[index].number_fluid < neighbors->max_neighbors) {
                neighbors->hash[index].fluid_particles[neighbors->hash[index].number_fluid] = p;
                neighbors->hash[index].number_fluid++;
            }
        }

        // Second pass - fill particle neighbors by processing buckets
      	// Could also iterate through hash directly but particles array will be shorter
        for (i=0; i<n_f; i++) {
           // Calculate hash index of bucket
           p = &fluid_particles[i];
           px = p->x_star;
           py = p->y_star;
           pz = p->z_star;

           index = hash_val(neighbors, px, py, pz);

           // If this bucket has been done try the next one
           if(neighbors->hash[index].hashed)
               continue;

            // Check neighbors of current bucket
            for (dx=-1; dx<=1; dx++) {
                x = px + dx*spacing;
                for (dy=-1; dy<=1; dy++) {
                    y = py + dy*spacing;
                    for (dz=-1; dz<=1; dz++) {
                        z = pz + dz*spacing;

                        // Make sure that the position is valid
                        if( floor(x/spacing) > neighbors->hash_size_x-1 || x < 0 ||
                            floor(y/spacing) > neighbors->hash_size_y-1 || y < 0 ||
                            floor(z/spacing) > neighbors->hash_size_z-1 || z < 0)
                          continue;

                        // Calculate hash index at neighbor point
                        neighbor_index = hash_val(neighbors, x, y, z);

                        // Add neighbor particles to each particle in current bucket
                        for (c=0; c<neighbors->hash[index].number_fluid; c++) {
			                // Particle in currently being worked on bucket
                            q = neighbors->hash[index].fluid_particles[c];
                            ne = &neighbors->particle_neighbors[q->id];
			                for(n=0; n<neighbors->hash[neighbor_index].number_fluid; n++){
                            // Append neighbor to q's neighbor list
		   	                    q_neighbor = neighbors->hash[neighbor_index].fluid_particles[n];
                                if(q->id == q_neighbor->id)
                                     continue;

                                r = sqrt((q_neighbor->x_star-q->x_star)*(q_neighbor->x_star-q->x_star)
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
