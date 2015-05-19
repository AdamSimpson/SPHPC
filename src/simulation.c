#include "simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "debug.h"
#include "communication.h"
#include "neighbors.h"
#include "fluid.h"
#include "geometry.h"
#include "fileio.h"
#include "input_parser.h"
#include "mpi.h"

int main(int argc, char *argv[]) {
  InitCommunication(argc, argv);

  struct Params params;
  struct AABB water_volume_global;
  struct AABB boundary_global;
  struct Communication communication;
  struct FluidParticles particles;
  struct Neighbors neighbors;
  struct FileIO file_io;

  SetParameters(&params, &neighbors, &boundary_global, &water_volume_global);

  SetParticleNumbers(&water_volume_global, &params, &communication);

  AllocateFluid(&particles, &params);

  AllocateNeighbors(&neighbors, &params, &boundary_global);

  AllocateCommunication(&communication);

  FileIOInit(&file_io, &params);

  // Initialize particles
  InitParticles(&particles, &params, &water_volume_global);

  // Print some parameters
  printf("Rank: %d, particles: %d, smoothing radius: %f \n",
         params.rank, params.number_particles_local, params.smoothing_radius);

  // Initial configuration
  WriteMPI(&particles, &params, &file_io);

  MPI_Barrier(MPI_COMM_WORLD);
  const double start_time = MPI_Wtime();

  for (int n=0; n<params.number_steps; n++) {

    DEBUG_PRINT("Rank %d Entering fluid step %d with %d particles\n",
                params.rank, n, params.number_particles_local);

    ApplyGravity(&particles, &params);

    // Advance to predicted position
    PredictPositions(&particles, &params, &boundary_global);

    if (n % 10 == 0)
      CheckPartition(&params);

    // Identify out of bounds particles and send them to appropriate rank
    IdentifyOOBParticles(&particles, &params, &communication);

    HaloExchange(&communication, &params, &particles);

    FindAllNeighbors(&params, &particles, &neighbors);

    const int solve_iterations = 4;
    int si;
    for (si=0; si<solve_iterations; si++) {
      ComputeDensities(&particles, &params, &neighbors);

      CalculateLambda(&particles, &params, &neighbors);

      UpdateHaloLambdas(&communication, &params, &particles);

      UpdateDPs(&particles, &params, &neighbors);

      UpdatePositionStars(&particles, &params, &boundary_global);

      UpdateHaloPositions(&communication, &params, &particles);
    }

    UpdateVelocities(&particles, &params);

    XSPHViscosity(&particles, &params, &neighbors);

    VorticityConfinement(&particles, &params, &neighbors);

    UpdatePositions(&particles, &params);

    // Write file at 30 FPS
    if (n % (int)(1.0/(params.time_step*30.0)) )
      WriteMPI(&particles, &params, &file_io);

  }

  const double end_time = MPI_Wtime();
  printf("Rank %d Elapsed seconds: %f, num particles: %d\n",
         params.rank, end_time-start_time, params.number_particles_local);

  // Release memory
  FreeFluid(&particles);
  FreeCommunication(&communication);
  FreeNeighbors(&neighbors);

  FileIOFinalize(&file_io);

  // Close MPI
  FinalizeCommunication();

  return 0;
}

void SetParameters(struct Params *const params,
                   struct Neighbors *const neighbors,
                   struct AABB *const boundary_global,
                   struct AABB *const water_volume_global)
{
  params->rank = get_rank();
  params->num_procs = get_num_procs();

  ReadParameters(params, boundary_global, water_volume_global);

  // Cubed volume
  const double volume = (water_volume_global->max_x - water_volume_global->min_x)
                      * (water_volume_global->max_y - water_volume_global->min_y)
                      * (water_volume_global->max_z - water_volume_global->min_z);

  // Initial spacing between particles
  const float spacing_particle = pow(volume/params->number_particles_global,
                                     1.0/3.0);

  // Let mass of each particle equal 1
  params->rest_density = params->number_particles_global/volume;
  printf("rest density: %f\n", params->rest_density);

  // Smoothing radius, h
  params->smoothing_radius = 2.0*spacing_particle;
  params->dq = 0.1*params->smoothing_radius;

  printf("smoothing radius: %f\n", params->smoothing_radius);

  // Set initial node boundaries
  // Equally divide water volume between all ranks
  double water_length = water_volume_global->max_x - water_volume_global->min_x;
  double equal_spacing =  water_length / params->num_procs;
  params->node_start_x = water_volume_global->min_x
                       + params->rank * equal_spacing;

  params->node_end_x = params->node_start_x + equal_spacing;

  // "Stretch" first and last ranks bounds to match boundary
  if (params->rank == 0)
    params->node_start_x  = boundary_global->min_x;
  if (params->rank == params->num_procs-1)
    params->node_end_x   = boundary_global->max_x;

  neighbors->max_neighbors = 60;
  neighbors->hash_spacing = params->smoothing_radius;
}
