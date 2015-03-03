#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "communication.h"
#include "fluid.h"
#include "geometry.h"
#include "fileio.h"
#include "simulation.h"

int main(int argc, char *argv[])
{
    init_communication(argc, argv);

    param_t params;
    AABB_t water_volume_global;
    AABB_t boundary_global;
    edge_t edges;
    oob_t out_of_bounds;
    fluid_particle_t *fluid_particles = NULL;
    neighbors_t neighbors;

    set_parameters(&params, &neighbors, &boundary_global, &water_volume_global);

    setParticleNumbers(&boundary_global, &water_volume_global, &edges, &out_of_bounds, &params);

    allocate_fluid(&fluid_particles, &params);

    allocate_neighbors(&neighbors, &boundary_global, &params);

    allocate_communication(&edges, &out_of_bounds);

    // Initialize particles
    initParticles(fluid_particles, &neighbors, &water_volume_global, &boundary_global, &edges, &params);

    // Print some parameters
    printf("Rank: %d, fluid_particles: %d, smoothing radius: %f \n", params.rank, params.number_fluid_particles_local, params.smoothing_radius);

    // Initial configuration
    int fileNum=0;
    writeMPI(fluid_particles, fileNum++, &params);

    // Main loop
    int n;
    double start_time, end_time;

    // Number of steps before frame needs to be written for 30 fps
    int steps_per_frame = (int)(1.0/(params.time_step*30.0));

    MPI_Barrier(MPI_COMM_WORLD);

    start_time = MPI_Wtime();
    for(n=0; n<params.number_steps; n++) {

        printf("Rank %d Entering fluid step %d with %d particles\n",params.rank, n, params.number_fluid_particles_local);

        apply_gravity(fluid_particles, &params);

        // Advance to predicted position
        predict_positions(fluid_particles, &boundary_global, &params);

        // Check boundary partitions
        if (n % 10 == 0)
            checkPartition(fluid_particles, &out_of_bounds, &params);

        // Identify out of bounds particles and send them to appropriate rank
        identify_oob_particles(fluid_particles, &out_of_bounds, &params);

        startHaloExchange(fluid_particles, &edges, &params);

        hash_fluid(fluid_particles, &neighbors, &boundary_global, &params);

        finishHaloExchange(fluid_particles, &edges, &params);

        hash_halo(fluid_particles, &neighbors, &boundary_global, &params);

        int solve_iterations = 4;
        int si;
        for(si=0; si<solve_iterations; si++)
        {
            compute_densities(fluid_particles, &neighbors, &params);

            calculate_lambda(fluid_particles, &neighbors, &params);

            update_halo_lambdas(fluid_particles, &edges, &params);

            update_dp(fluid_particles, &neighbors, &params);

            update_dp_positions(fluid_particles, &boundary_global, &params);

            update_halo_positions(fluid_particles, &edges, &params);
        }

        update_velocities(fluid_particles, &params);

        XSPH_viscosity(fluid_particles, &neighbors, &params);

        vorticity_confinement(fluid_particles, &neighbors, &params);

        update_positions(fluid_particles, &params);

        if (n % steps_per_frame == 0)
            writeMPI(fluid_particles, fileNum++, &params);

    }
    end_time = MPI_Wtime();
    printf("Rank %d Elapsed seconds: %f, num particles: %d\n", params.rank, end_time-start_time, params.number_fluid_particles_local);

    // Release memory
    // Need to fix this to free everything correctly
    free(fluid_particles);
    free(edges.edge_indices_left);
    free(edges.edge_indices_right);
    free(out_of_bounds.oob_indices_left);
    free(out_of_bounds.oob_indices_right);

    // Close MPI
    freeMpiTypes();
    MPI_Finalize();

    return 0;
}

void set_parameters(param_t *params, neighbors_t *neighbors,
                    AABB_t *boundary_global, AABB_t* water_volume_global)
{
  params->rank = get_rank();
  params->nprocs = get_num_procs();
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
  double volume = (water_volume_global->max_x - water_volume_global->min_x)
                * (water_volume_global->max_y - water_volume_global->min_y)
                * (water_volume_global->max_z - water_volume_global->min_z);

  // Initial spacing between particles
  float spacing_particle = pow(volume/params->number_fluid_particles_global,1.0/3.0);

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
  double equal_spacing =  water_length / params->nprocs;
  params->node_start_x = water_volume_global->min_x
                      + params->rank * equal_spacing;

  params->node_end_x = params->node_start_x + equal_spacing;

  // "Stretch" first and last ranks bounds to match boundary
  if (params->rank == 0)
      params->node_start_x  = boundary_global->min_x;
  if (params->rank == params->nprocs-1)
      params->node_end_x   = boundary_global->max_x;

    neighbors->max_neighbors = 60;
    neighbors->hash_spacing = params->smoothing_radius;
}
