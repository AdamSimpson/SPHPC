#include "simulation.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "debug.h"
#include "communication.h"
#include "neighbors.h"
#include "particles.h"
#include "geometry.h"
#include "fileio.h"
#include "input_parser.h"

int main(int argc, char *argv[]) {
  InitCommunication(argc, argv);

  struct Params params;
  struct AABB water_volume_global;
  struct AABB boundary_global;
  struct Communication communication;
  struct Particles particles;
  struct Neighbors neighbors;
  struct FileIO file_io;

  SetParameters(&params, &particles, &neighbors,
                &boundary_global, &water_volume_global);

  SetParticleNumbers(&water_volume_global, &particles, &params, &communication);

  AllocInitParticles(&particles, &params, &water_volume_global);

  AllocInitNeighbors(&neighbors, &particles, &params, &boundary_global);

  AllocateCommunication(&communication);

  AllocInitFileIO(&file_io, &particles, &params);

  WriteParticles(&particles, &params, &file_io);

  // Copy over structs with only scalar values
  #pragma acc enter data copyin(params, water_volume_global, boundary_global)

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

      UpdatePositionStars(&particles, &boundary_global);
      UpdateHaloTuple(&communication, &params, &particles,
                       particles.x_star, particles.y_star, particles.z_star);
    }

    UpdateVelocities(&particles, &params);

//    ApplyViscosity(&particles, &params, &neighbors);

    ApplyVorticityConfinement(&particles, &params, &neighbors);

    UpdatePositions(&particles);

    // Write file at 30 FPS and print average density
    if (n % (int)(1.0/(params.time_step*30.0)) ) {
      #pragma acc update host(particles.x[:particles.local_count], particles.y[:particles.local_count], particles.z[:particles.local_count])
      WriteParticles(&particles, &params, &file_io);
//      PrintAverageDensity(&particles);
    }
  }

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
                   struct AABB *const water_volume_global) {
  params->rank = get_rank();
  params->proc_count = get_proc_count();

  ReadParameters(params, particles, boundary_global, water_volume_global);

  // Cubed volume
  const double volume = (water_volume_global->max_x - water_volume_global->min_x)
                      * (water_volume_global->max_y - water_volume_global->min_y)
                      * (water_volume_global->max_z - water_volume_global->min_z);

  // Initial spacing between particles
  const float spacing_particle = pow(volume/particles->global_count, 1.0/3.0);

  // Let mass of each particle equal 1
  params->rest_density = particles->global_count/volume;
  printf("rest density: %f\n", params->rest_density);

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

  params->W_norm = 315.0/(64.0*M_PI*pow(params->smoothing_radius, 9.0));
  params->DelW_norm = -45.0/(M_PI*pow(params->smoothing_radius, 6.0));
}
