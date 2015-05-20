////////////////////////////////////////////////
// Utility Functions
////////////////////////////////////////////////
#include "fileio.h"
#include "simulation.h"
#include "debug.h"
#include "safe_alloc.h"
#include "fluid.h"
#include "mpi.h"
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

void FileIOInit(struct FileIO *const file_io,
                const struct Params *const params) {
  int num_components = 3;
  file_io->write_buffer = SAFE_ALLOC(num_components*params->max_particles_local,
                                     sizeof(double));

  const char output_path[] = "/home/atj/OLCF/SPHPC/output";
  file_io->output_path = SAFE_ALLOC(strlen(output_path)+2, sizeof(char));
  strcpy(file_io->output_path, output_path);

  file_io->file_num = 0;
}

void FileIOFinalize(struct FileIO *const file_io) {
    free(file_io->write_buffer);
    free(file_io->output_path);
}

// Write boundary in MPI
void WriteMPI(const struct FluidParticles *const particles,
              const struct Params *const params,
              struct FileIO *const file_io) {

  MPI_File file;
  MPI_Status status;
  char file_name[128];
  sprintf(file_name, "%s/sim-%d.bin", file_io->output_path, file_io->file_num);

  const int num_particles = params->number_particles_local;

  // How many bytes each process will write
  int rank_write_counts[params->num_procs];
  // alltoall of write counts
  int num_doubles_to_send = 3 * num_particles;
  MPI_Allgather(&num_doubles_to_send, 1, MPI_INT, rank_write_counts,
                1, MPI_INT, MPI_COMM_WORLD);

  // Displacement can overflow with int, max size = 8*3*(global num particles)
  uint64_t displacement=0;
  for (int i=0; i<params->rank; i++)
    displacement+=rank_write_counts[i];
  // Displacement is in bytes
  displacement*=sizeof(double);

  int index=0;
  for (int i=0; i<num_particles; i++) {
    file_io->write_buffer[index]   = particles->x[i];
    file_io->write_buffer[index+1] = particles->y[i];
    file_io->write_buffer[index+2] = particles->z[i];
    index+=3;
  }

  // Open file
  MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &file);

  // write coordinates to file
  MPI_File_write_at(file, displacement, file_io->write_buffer,
                    num_doubles_to_send, MPI_DOUBLE, &status);

  int num_bytes;
  MPI_Get_elements(&status, MPI_CHAR, &num_bytes);
  DEBUG_PRINT("rank %d wrote %d bytes(%d particles) to file %s\n",
              params->rank,num_bytes,num_particles,name);

  // Close file
  MPI_File_close(&file);

  ++file_io->file_num;
}
