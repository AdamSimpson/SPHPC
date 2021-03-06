////////////////////////////////////////////////
// Utility Functions
////////////////////////////////////////////////
#include "fileio_MPI_IO.h"
#include "simulation.h"
#include "print_macros.h"
#include "safe_alloc.h"
#include "particles.h"
#include "mpi.h"
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
  int num_components = 3;
  file_io->write_buffer = SAFE_ALLOC(num_components*particles->max_local,
                                     sizeof(double));

  #pragma acc enter data copyin( file_io[:1],                    \
    file_io->write_buffer[0:num_components*particles->max_local])

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
}

void FinalizeFileIO(struct FileIO *const file_io) {
//  #pragma acc exit data delete(file_io->write_buffer)
  free(file_io->write_buffer);
  free(file_io->output_path);
}

// Write particle coordinates in MPI
void WriteParticles(const struct Particles *const particles,
              const struct Params *const params,
              struct FileIO *const file_io) {

  MPI_File file;
  MPI_Status status;
  char file_name[1024];
  sprintf(file_name, "%s/sim-%d.bin", file_io->output_path, file_io->file_num);

  const int num_particles = particles->local_count;

  double *x = particles->x;
  double *y = particles->y;
  double *z = particles->z;
  double *write_buffer = file_io->write_buffer;

  #pragma acc parallel loop vector_length(1024) present(x,y,z,write_buffer)
  for (int i=0; i<num_particles; ++i) {
    write_buffer[3*i]   = x[i];
    write_buffer[3*i+1] = y[i];
    write_buffer[3*i+2] = z[i];
  }
  #pragma acc update host(write_buffer[:num_particles*3])

  // How many bytes each process will write
  int rank_write_counts[params->proc_count];
  // alltoall of write counts
  int num_doubles_to_send = 3 * num_particles;
  MPI_Allgather(&num_doubles_to_send, 1, MPI_INT, rank_write_counts,
                1, MPI_INT, MPI_COMM_WORLD);

  // Displacement can overflow with int, max size = 8*3*(global num particles)
  uint64_t displacement=0;
  for (int i=0; i<params->rank; ++i)
    displacement+=rank_write_counts[i];
  // Displacement is in bytes
  displacement*=sizeof(double);

  // Open file
  MPI_File_open(MPI_COMM_WORLD, file_name, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &file);

  // write coordinates to file
  MPI_File_write_at(file, displacement, file_io->write_buffer,
                    num_doubles_to_send, MPI_DOUBLE, &status);

  int num_bytes;
  MPI_Get_elements(&status, MPI_CHAR, &num_bytes);
  DEBUG_PRINT("rank %d wrote %d bytes(%d particles) to output\n",
              params->rank,num_bytes,num_particles);

  // Close file
  MPI_File_close(&file);

  ++file_io->file_num;
}
