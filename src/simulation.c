#include "simulation.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "communication.h"
#include "fluid.h"
#include "geometry.h"
#include "fileio.h"

int main(int argc, char *argv[])
{
  InitCommunication(argc, argv);

  Params params;
  AABB water_volume_global;
  AABB boundary_global;
  Communication communication;
  FluidParticle *fluid_particles = NULL;
  Neighbors neighbors;

  SetParameters(&params, &neighbors, &boundary_global, &water_volume_global);

  SetParticleNumbers(&water_volume_global, &params, &communication);

  AllocateFluid(&fluid_particles, &params);

  AllocateNeighbors(&neighbors, &params, &boundary_global);

  AllocateCommunication(&communication);

  // Initialize particles
  InitParticles(fluid_particles, &params, &water_volume_global);

  // Print some parameters
  printf("Rank: %d, fluid_particles: %d, smoothing radius: %f \n", params.rank, params.number_fluid_particles_local, params.smoothing_radius);

  // Initial configuration
  int fileNum=0;
  WriteMPI(fluid_particles, &params, fileNum++);

  MPI_Barrier(MPI_COMM_WORLD);
  const double start_time = MPI_Wtime();

  for (int n=0; n<params.number_steps; n++) {

    printf("Rank %d Entering fluid step %d with %d particles\n",params.rank, n, params.number_fluid_particles_local);

    ApplyGravity(fluid_particles, &params);

    // Advance to predicted position
    PredictPositions(fluid_particles, &params, &boundary_global);

    if (n % 10 == 0)
      CheckPartition(&params);

    // Identify out of bounds particles and send them to appropriate rank
    IdentifyOOBParticles(fluid_particles, &params, &communication);

    StartHaloExchange(&communication, &params, fluid_particles);

    HashFluid(fluid_particles, &params, &boundary_global, &neighbors);

    FinishHaloExchange(&communication, fluid_particles, &params);

    HashHalo(fluid_particles, &params, &boundary_global, &neighbors);

    const int solve_iterations = 4;
    int si;
    for (si=0; si<solve_iterations; si++) {
      ComputeDensities(fluid_particles, &params, &neighbors);

      CalculateLambda(fluid_particles, &params, &neighbors);

      UpdateHaloLambdas(&communication, &params, fluid_particles);

      UpdateDPs(fluid_particles, &params, &neighbors);

      UpdatePositionStars(fluid_particles, &params, &boundary_global);

      UpdateHaloPositions(&communication, &params, fluid_particles);
    }

    UpdateVelocities(fluid_particles, &params);

    XSPHViscosity(fluid_particles, &params, &neighbors);

    VorticityConfinement(fluid_particles, &params, &neighbors);

    UpdatePositions(fluid_particles, &params);

    // Write file at 30 FPS
    if (n % (int)(1.0/(params.time_step*30.0)) )
      WriteMPI(fluid_particles, &params, fileNum++);

  }

  const double end_time = MPI_Wtime();
  printf("Rank %d Elapsed seconds: %f, num particles: %d\n", params.rank, end_time-start_time, params.number_fluid_particles_local);

  // Release memory
  free(fluid_particles);
  FreeCommunication(&communication);
  FreeNeighbors(&neighbors);

  // Close MPI
  FinalizeCommunication();

  return 0;
}

void SetParameters(Params *const params,
                   Neighbors *const neighbors,
                   AABB *const boundary_global,
                   AABB *const water_volume_global)
{
  params->rank = get_rank();
  params->num_procs = get_num_procs();
  params->g = 9.8;
  params->number_steps = 2000;
  params->time_step = 1.0/60.0;
  params->c = 0.01;
  params->k = 0.05;
  // Approximate
  params->number_fluid_particles_global = 65536*2;

  boundary_global->min_x = 0.0;
  boundary_global->max_x = 100.0;
  boundary_global->min_y = 0.0;
  boundary_global->max_y = 80.0;
  boundary_global->min_z = 0.0;
  boundary_global->max_z = 30.0;

  water_volume_global->min_x = 0.1;
  water_volume_global->max_x = boundary_global->max_x - 20.0;
  water_volume_global->min_y = 0.1;
  water_volume_global->max_y = boundary_global->max_y - 30.0;
  water_volume_global->min_z = 0.1;
  water_volume_global->max_z = boundary_global->max_z - 10.0;

  // Cubed volume
  const double volume = (water_volume_global->max_x - water_volume_global->min_x)
                      * (water_volume_global->max_y - water_volume_global->min_y)
                      * (water_volume_global->max_z - water_volume_global->min_z);

  // Initial spacing between particles
  const float spacing_particle = pow(volume/params->number_fluid_particles_global,1.0/3.0);

  // Let mass of each particle equal 1
  params->rest_density = params->number_fluid_particles_global/volume;
  printf("rest density: %f\n", params->rest_density);

  // Smoothing radius, h
  params->smoothing_radius = 2.0*spacing_particle;
  params->dq = 0.1*params->smoothing_radius;

  printf("smoothing radius: %f\n", params->smoothing_radius);

  // Set initial node boundaries
  // Equally divide water volume between all ranks
  double water_length = (water_volume_global->max_x - water_volume_global->min_x);
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
