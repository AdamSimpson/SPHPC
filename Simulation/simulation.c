#include "simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "print_macros.h"
#include "communication.h"
#include "neighbors.h"
#include "particles.h"
#include "geometry.h"
#include "obstacle.h"
#include "input_parser.h"
#ifdef USE_ADIOS
#include "fileio_adios.h"
#else
#include "fileio_MPI_IO.h"
#endif


int main(int argc, char *argv[]) {
  InitCommunication(argc, argv);

  struct Params params;
  struct AABB water_volume_global;
  struct AABB boundary_global;
  struct Communication communication;
  struct Particles particles;
  struct Neighbors neighbors;
  struct FileIO file_io;
  struct Obstacle obstacle;
  struct Obstacle_CUDA obstacle_cuda;

  SetParameters(&params, &particles, &neighbors,
                &boundary_global, &water_volume_global, &obstacle);

  SetParticleNumbers(&water_volume_global, &particles, &params, &communication);

  AllocInitParticles(&particles, &params, &water_volume_global);

  AllocInitNeighbors(&neighbors, &particles, &params, &boundary_global);

  AllocateCommunication(&communication);

  AllocInitObstacle(&obstacle);

  AllocInitFileIO(&file_io, &particles, &params);

  WriteParticles(&particles, &params, &file_io);

  // Copy over structs with only scalar values
  #pragma acc enter data copyin(params, water_volume_global, boundary_global)

  double start_time = GetTime();

  for (int n=0; n<params.number_steps; ++n) {

    DEBUG_PRINT("Rank %d Entering fluid step %d with %d particles\n",
                params.rank, n, particles.local_count);

    ApplyGravity(&particles, &params);

    PredictPositions(&particles, &params, &boundary_global);

    if (n % 10 == 0)
      BalanceNodes(&particles, &params);

    ExchangeOOB(&communication, &particles, &params);

    ExchangeHalo(&communication, &params, &particles);

    FindAllNeighbors(&particles, &params, &neighbors);

    const int solve_iterations = 4;
    for (int sub_i=0; sub_i<solve_iterations; ++sub_i) {
      ComputeDensities(&particles, &params, &neighbors);
      UpdateHaloScalar(&communication, &params, &particles, particles.density);

      ComputeLambda(&particles, &params, &neighbors);
      UpdateHaloScalar(&communication, &params, &particles, particles.lambda);

      UpdateDPs(&particles, &params, &neighbors);
      UpdateHaloTuple(&communication, &params, &particles,
                       particles.dp_x, particles.dp_y, particles.dp_z);

      UpdatePositionStars(&particles, &boundary_global, &obstacle, &params);
      UpdateHaloTuple(&communication, &params, &particles,
                       particles.x_star, particles.y_star, particles.z_star);

      // Put collision detection here explicity
    }

    UpdateVelocities(&particles, &params);

    ApplyViscosity(&particles, &params, &neighbors);
    UpdateHaloTuple(&communication, &params, &particles,
                    particles.v_x, particles.v_y, particles.v_z);

    ComputeVorticity(&particles, &params, &neighbors);
    UpdateHaloTuple(&communication, &params, &particles,
                    particles.w_x, particles.w_y, particles.w_z);
    ApplyVorticityConfinement(&particles, &params, &neighbors);

    UpdatePositions(&particles);

    PrintAverageDensity(&particles);

    // Write file at 30 FPS
    if (n % (int)(1.0/(params.time_step*30.0)) == 0)
      WriteParticles(&particles, &params, &file_io);

  }

  printf("Total loop time: %f seconds\n", GetTime() - start_time);

  FinalizeParticles(&particles);
  FinalizeNeighbors(&neighbors);
  FinalizeFileIO(&file_io);

  FreeCommunication(&communication);
  FinalizeCommunication();

  return 0;
}

void SetParameters(struct Params *const params,
                   struct Particles *const particles,
                   struct Neighbors *const neighbors,
                   struct AABB *const boundary_global,
                   struct AABB *const water_volume_global,
                   struct Obstacle *const obstacle) {
  params->rank = get_rank();
  params->proc_count = get_proc_count();

  ReadParameters(params, particles, boundary_global, water_volume_global, obstacle);

  // Initial liquid volume
  const double volume = (water_volume_global->max_x - water_volume_global->min_x)
                      * (water_volume_global->max_y - water_volume_global->min_y)
                      * (water_volume_global->max_z - water_volume_global->min_z);

  // Initial spacing between particles
  const float spacing_particle = pow(volume/particles->global_count, 1.0/3.0);

  // Smoothing radius, h
  params->smoothing_radius = 2.0*spacing_particle;
  params->dq = 0.1*params->smoothing_radius;

  printf("smoothing radius: %f\n", params->smoothing_radius);

  // Set initial node boundaries
  // Equally divide water volume between all ranks
  double water_length = water_volume_global->max_x - water_volume_global->min_x;
  double equal_spacing =  water_length / params->proc_count;
  params->node_start_x = water_volume_global->min_x
                       + params->rank * equal_spacing;

  params->node_end_x = params->node_start_x + equal_spacing;

  // "Stretch" first and last ranks bounds to match boundary
  if (params->rank == 0)
    params->node_start_x  = boundary_global->min_x;
  if (params->rank == params->proc_count-1)
    params->node_end_x   = boundary_global->max_x;

  #ifndef M_PI
    #define M_PI  3.14159265358979323846
  #endif

  // Normalization for poly6 and del_spiky uing W(q) where q = r/h
  params->W_norm = 315.0/(64.0*M_PI*pow(params->smoothing_radius, 3.0));
  params->DelW_norm = -45.0/(M_PI*pow(params->smoothing_radius, 3.0));
}
