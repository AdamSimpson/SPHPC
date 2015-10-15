////////////////////////////////////////////////
// Utility Functions
////////////////////////////////////////////////
#include "fileio_adios.h"
#include "simulation.h"
#include "print_macros.h"
#include "safe_alloc.h"
#include "particles.h"
#include "mpi.h"
#include "adios.h"
#include <stdio.h>
#include <errno.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

void AllocInitFileIO(struct FileIO *const file_io,
                const struct Particles *const particles,
                const struct Params *params) {
  // Get current path
  char current_path[1024];
  getcwd(current_path, sizeof(current_path));
  if(current_path == NULL) {
    printf("Error getting current working directory \n");
    exit(-1);
  }

  char output_path[1024];
  sprintf(output_path, "%s/output", current_path);

  // Create pwd/output directory
  if (params->rank == 0) {
    int dir_return = mkdir(output_path, 0775);
//    if (dir_return != 0) {
//      printf("Error creating output directory: %s\n", strerror(errno));
//      exit(-1);
//    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Set output directory
  file_io->output_path = SAFE_ALLOC(strlen(output_path)+2, sizeof(char));
  strcpy(file_io->output_path, output_path);

  file_io->file_num = 0;

  adios_init ("sph_config.xml", MPI_COMM_WORLD);

  file_io->adios_handle = NULL;

  file_io->rank = params->rank;
}

void FinalizeFileIO(struct FileIO *const file_io) {
  free(file_io->output_path);
  adios_finalize(file_io->rank);
}

// Write particle coordinates with adios
void WriteParticles(const struct Particles *const particles,
              const struct Params *const params,
              struct FileIO *const file_io) {

  char file_name[1024];
  sprintf(file_name, "%s/sim-%d.bp", file_io->output_path, file_io->file_num);

  int local_count = particles->local_count;
  int global_count = particles->global_count;

  // Update x,y,z components on host
  #pragma acc update host(particles->x[:local_count],  \
                          particles->y[:local_count],  \
                          particles->z[:local_count],  \
                          particles->density[:local_count]) async

  // Write positions to file
  int adios_err;
  uint64_t adios_groupsize, adios_totalsize;

  adios_open(&file_io->adios_handle, "particles", file_name, "w", MPI_COMM_WORLD);

  // Compute adios local size to write
  adios_groupsize = 3 * sizeof(int)
                  + 1 * sizeof(double)
                  + 4 * local_count * sizeof(double);  // (x,y,z,density)
  adios_group_size(file_io->adios_handle, adios_groupsize, &adios_totalsize);

  // Compute offset in global output for current rank(sum of ranks to the "left")
  int offset_count = 0;
  MPI_Exscan(&local_count, &offset_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // Wait for host data update to complete
  #pragma acc wait

  adios_write(file_io->adios_handle, "global_particle_count", &global_count);
  adios_write(file_io->adios_handle, "local_particle_count", &local_count);
  adios_write(file_io->adios_handle, "offset_particle_count", &offset_count);
  adios_write(file_io->adios_handle, "particle_radius", &(particles->rest_radius));
  adios_write(file_io->adios_handle, "x", particles->x);
  adios_write(file_io->adios_handle, "y", particles->y);
  adios_write(file_io->adios_handle, "z", particles->z);
  adios_write(file_io->adios_handle, "density", particles->density);

  adios_close(file_io->adios_handle);

  ++file_io->file_num;
 }

// Write color field with adios
void WriteColorField(const struct Particles *const particles,
                     const struct Params *const params,
                     struct FileIO *const file_io) {
/*
  char file_name[1024];
  sprintf(file_name, "%s/sim-%d.bp", file_io->output_path, file_io->file_num);

  int particle_count = particles->local_count;

  // Update x,y,z components on host
  #pragma acc update host(particles->x[:particle_count],  \
                          particles->y[:particle_count],  \
                          particles->z[:particle_count]) async

  // Write positions to file
  int adios_err;
  uint64_t adios_groupsize, adios_totalsize;

  adios_open(&file_io->adios_handle, "color_field", file_name, "w", MPI_COMM_WORLD);

  adios_groupsize = sizeof(int)                            // component_count
                   + 3 * particle_count * sizeof(double);  // (x,y,z) coords
  adios_group_size(file_io->adios_handle, adios_groupsize, &adios_totalsize);

  // Wait for host data update to complete
  #pragma acc wait

  adios_write(file_io->adios_handle, "particle_count", &particle_count);
  adios_write(file_io->adios_handle, "x", particles->x);
  adios_write(file_io->adios_handle, "y", particles->y);
  adios_write(file_io->adios_handle, "z", particles->z);

  adios_close(file_io->adios_handle);

  ++file_io->file_num;
  */
}
