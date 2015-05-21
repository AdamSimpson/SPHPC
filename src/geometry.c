#include "geometry.h"
#include <stdio.h>
#include <math.h>
#include "particles.h"
#include "communication.h"
#include "simulation.h"
#include "debug.h"
#include "mpi.h"

////////////////////////////////////////////////
// Utility Functions
////////////////////////////////////////////////
static inline double min(const double a, const double b) {
  double min = a;
  min = b < min ? b : min;
  return min;
}

void ConstructFluidVolume(struct Particles *const particles,
                          struct Params *const params,
                          const struct AABB *const fluid) {

  const double spacing = params->smoothing_radius/2.0;

  // Start node particles at integer multiple of spacing
  double min_x;
  if (params->rank == 0)
    min_x = fluid->min_x;
  else {
    int num_to_left = floor((params->node_start_x - fluid->min_x)/spacing);
    min_x = num_to_left*spacing + fluid->min_x;
  }

  const double max_x = min(fluid->max_x, params->node_end_x);

  const int num_x = floor((max_x - min_x) / spacing);
  const int num_y = floor((fluid->max_y - fluid->min_y ) / spacing);
  const int num_z = floor((fluid->max_z - fluid->min_z ) / spacing);

  // Place particles inside bounding volume
  int i = 0;
  for (int nz=0; nz<num_z; nz++) {
    const double z = fluid->min_z + nz*spacing + spacing/2.0;
    for (int ny=0; ny<num_y; ny++) {
      const double y = fluid->min_y + ny*spacing + spacing/2.0;
      for(int nx=0; nx<num_x; nx++) {
        const double x = min_x + nx*spacing + spacing/2.0;
        particles->x[i] = x;
        particles->y[i] = y;
        particles->z[i] = z;
        particles->id[i] = i;
        i++;
      }
    }
  }

  particles->local_count = i;
  DEBUG_PRINT("rank %d: min fluid: %f max fluid x: %f\n",
              params->rank, min_x + spacing/2.0,
              min_x + num_x*spacing + spacing/2.0);
}

// Sets upper bound on number of particles, used for memory allocation
void SetParticleNumbers(const struct AABB *const fluid_global,
                        struct Particles *particles,
                        struct Params *const params,
                        struct Communication *const communication) {

  const double spacing = params->smoothing_radius/2.0;

  // Get some baseline numbers useful to define maximum particle numbers
  const int num_x = floor((fluid_global->max_x - fluid_global->min_x ) / spacing);
  const int num_y = floor((fluid_global->max_y - fluid_global->min_y ) / spacing);
  const int num_z = floor((fluid_global->max_z - fluid_global->min_z ) / spacing);

  // Initial fluid particles per rank
  const int num_initial = (num_x * num_y * num_z)/get_proc_count();

  // Set max communication particles to a 10th of node start number
  communication->max_particles = num_initial/10;

  // Add initial space and left/right out of bounds/halo particles
  particles->max_local = num_initial + 4*communication->max_particles;

  DEBUG_PRINT("initial number of particles %d\n", num_initial);
  DEBUG_PRINT("Max fluid particles local: %d\n", params->max_particles_local);
}

// Test if boundaries need to be adjusted
void CheckPartition(const struct Particles *particles,
                    struct Params *const params) {
  const int num_rank = particles->local_count;
  const int rank = params->rank;
  const int num_procs = params->proc_count;
  const double h = params->smoothing_radius;

  // Setup nodes to left and right of self
  const int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  const int proc_to_right = (rank == num_procs-1 ? MPI_PROC_NULL : rank+1);

  // Get number of particles and partition length  from right and left
  const double length = params->node_end_x - params->node_start_x;
  double node[2]  = {(double)num_rank, length};
  double left[2]  = {0.0, 0.0};
  double right[2] = {0.0, 0.0};
  int tag = 627;
  // Send number of particles to  right and receive from left
  MPI_Sendrecv(node, 2, MPI_DOUBLE, proc_to_right, tag,
               left, 2, MPI_DOUBLE, proc_to_left,  tag,
               MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  // Send number of particles to left and receive from right
  tag = 895;
  MPI_Sendrecv(node,  2, MPI_DOUBLE, proc_to_left, tag,
               right, 2, MPI_DOUBLE, proc_to_right, tag,
               MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  // Number of particles in left/right ranks
  const int num_right = (int)right[0];

  // Partition length of left/right ranks
  const double length_left = left[1];
  const double length_right = right[1];

  const int even_particles = particles->global_count
                            / (double)params->proc_count;
  const int max_diff = even_particles/10.0f;

  // Difference in particle numbers from an even distribution
  const int diff_right = num_right - even_particles;
  const int diff_self  = num_rank  - even_particles;

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
  if( diff_right > max_diff && length_right > 2*h && rank != num_procs-1)
    params->node_end_x += h;
  // Rank to right has too few particles and will move its start to right
  else if (diff_right < -max_diff && length > 2*h && rank != num_procs-1)
    params->node_end_x -= h;

  DEBUG_PRINT("rank %d node_start %f node_end %f \n",
              rank, params->node_start_x, params->node_end_x);
}
