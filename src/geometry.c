#include <stdio.h>
#include "geometry.h"
#include "communication.h"

void constructFluidVolume(fluid_particle_t *fluid_particles, AABB_t* fluid, param_t *params)
{
    double spacing, x, y, z, min_x, max_x;
    int num_x, num_y, num_z, nx, ny, nz;
    int i = 0;
    fluid_particle_t *p;

    spacing = params->smoothing_radius/2.0;

    // Start node particles at integer multiple of spacing
    if(params->rank == 0)
      min_x = fluid->min_x;
    else {
      int num_to_left = floor((params->node_start_x - fluid->min_x)/spacing);
      min_x = num_to_left*spacing + fluid->min_x;
    }

    max_x = min(fluid->max_x, params->node_end_x);

    num_x = floor((max_x - min_x) / spacing);
    num_y = floor((fluid->max_y - fluid->min_y ) / spacing);
    num_z = floor((fluid->max_z - fluid->min_z ) / spacing);

    // Place particles inside bounding volume
    for(nz=0; nz<num_z; nz++) {
        z = fluid->min_z + nz*spacing + spacing/2.0;
        for(ny=0; ny<num_y; ny++) {
            y = fluid->min_y + ny*spacing + spacing/2.0;
            for(nx=0; nx<num_x; nx++) {
                x = min_x + nx*spacing + spacing/2.0;
                p = &fluid_particles[i];
                p->x = x;
                p->y = y;
                p->z = z;
                p->id = i;
                i++;
            }
        }
    }
    params->number_fluid_particles_local = i;
    printf("rank %d: min fluid: %f max fluid x: %f\n", params->rank, min_x + spacing/2.0, x);
}

// Sets upper bound on number of particles, used for memory allocation
void setParticleNumbers(AABB_t *fluid_global, communication_t *communication, param_t *params)
{
    int num_x;
    int num_y;
    int num_z;
    double spacing = params->smoothing_radius/2.0;

    // Get some baseline numbers useful to define maximum particle numbers
    num_x = floor((fluid_global->max_x - fluid_global->min_x ) / spacing);
    num_y = floor((fluid_global->max_y - fluid_global->min_y ) / spacing);
    num_z = floor((fluid_global->max_z - fluid_global->min_z ) / spacing);

    // Initial fluid particles per rank
    int num_initial = (num_x * num_y * num_z)/get_num_procs();

    // Set max communication particles to a 10th of node start number
    communication->max_comm_particles = num_initial/10;

    // Add initial space and left/right out of bounds/halo particles
    params->max_fluid_particles_local = num_initial + 4*communication->max_comm_particles;

    printf("initial number of particles %d\n", num_initial);
    printf("Max fluid particles local: %d\n", params->max_fluid_particles_local);
}

// Test if boundaries need to be adjusted
void checkPartition(fluid_particle_t *fluid_particles, param_t *params)
{
    int num_rank = params->number_fluid_particles_local;
    int rank = params->rank;
    int nprocs = params->nprocs;
    double h = params->smoothing_radius;

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    // Get number of particles and partition length  from right and left
    double length = params->node_end_x - params->node_start_x;
    double node[2]  = {(double)num_rank, length};
    double left[2]  = {0.0, 0.0};
    double right[2] = {0.0, 0.0};
    int tag = 627;
    // Send number of particles to  right and receive from left
    MPI_Sendrecv(node, 2, MPI_DOUBLE, proc_to_right, tag, left,2,MPI_DOUBLE,proc_to_left,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    // Send number of particles to left and receive from right
    tag = 895;
    MPI_Sendrecv(node, 2, MPI_DOUBLE, proc_to_left, tag, right,2,MPI_DOUBLE,proc_to_right,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Number of particles in left/right ranks
    int num_left = (int)left[0];
    int num_right = (int)right[0];

    // Partition length of left/right ranks
    double length_left = left[1];
    double length_right = right[1];

    int even_particles = params->number_fluid_particles_global/(double)params->nprocs;
    int max_diff = even_particles/10.0f;

    // Difference in particle numbers from an even distribution
    int diff_left  = num_left  - even_particles;
    int diff_right = num_right - even_particles;
    int diff_self  = num_rank  - even_particles;

    // Look at "bins" formed by node start/ends from right to left
    // Only modify node start based upon bin to the left
    // Node end may be modified by node to the rights start
    // Particles per proc if evenly divided

    // current rank has too many particles
    if( diff_self > max_diff && length > 2*h && rank != 0)
        params->node_start_x += h;
    // current rank has too few particles
    else if (diff_self < -max_diff && length_left > 2*h && rank != 0)
        params->node_start_x -= h;

    // Rank to right has too many particles and with move its start to left
    if( diff_right > max_diff && length_right > 2*h && rank != nprocs-1)
        params->node_end_x += h;
    // Rank to right has too few particles and will move its start to right
    else if (diff_right < -max_diff && length > 2*h && rank != nprocs-1)
        params->node_end_x -= h;

    printf("rank %d node_start %f node_end %f \n", rank, params->node_start_x, params->node_end_x);
}

////////////////////////////////////////////////
// Utility Functions
////////////////////////////////////////////////
double min(double a, double b){
    double min = a;
    min = b < min ? b : min;
    return min;
}

double max(double a, double b){
    double max = a;
    max = b > max ? b : max;
    return max;
}

int sgn(double x) {
    int val = 0;
    if (x < 0.0)
        val = -1;
    else if (x > 0.0)
        val = 1;

    return val;
}
