////////////////////////////////////////////////
// Utility Functions
////////////////////////////////////////////////
#include "fileio.h"
#include "simulation.h"
#include "debug.h"
#include "safe_alloc.h"
#include "mpi.h"
#include <stdio.h>
#include <inttypes.h>
#include "stdlib.h"
#include "fluid.h"

// Write boundary in MPI
void WriteMPI(const struct FluidParticles *const particles,
              const struct Params *const params,
              const int fileNum) {

  MPI_File file;
  MPI_Status status;
  char name[64];
  sprintf(name, "/home/atj/OLCF/SPHPC/output/sim-%d.bin", fileNum);

  const int num_particles = params->number_particles_local;

  // How many bytes each process will write
  int rank_write_counts[params->num_procs];
  // alltoall of write counts
  int num_doubles_to_send = 3 * num_particles;
  MPI_Allgather(&num_doubles_to_send, 1, MPI_INT, rank_write_counts, 1, MPI_INT, MPI_COMM_WORLD);
  // Displacement can overflow with int, max size = 8*3*(global num particles)
  uint64_t displacement=0;
  for (int i=0; i<params->rank; i++)
    displacement+=rank_write_counts[i];
  // Displacement is in bytes
  displacement*=sizeof(double);

  // Write x,y,z data to memory
  double *const send_buffer = SAFE_ALLOC(num_doubles_to_send, sizeof(double));

  int index=0;
  for (int i=0; i<num_particles; i++) {
    send_buffer[index]   = particles->x[i];
    send_buffer[index+1] = particles->y[i];
    send_buffer[index+2] = particles->z[i];
    index+=3;
  }

  // Open file
  MPI_File_open(MPI_COMM_WORLD, name, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

  // write coordinates to file
  MPI_File_write_at(file, displacement, send_buffer , num_doubles_to_send, MPI_DOUBLE, &status);

  int num_bytes;
  MPI_Get_elements(&status, MPI_CHAR, &num_bytes);
  DEBUG_PRINT("rank %d wrote %d bytes(%d particles) to file %s\n",
              params->rank,num_bytes,num_particles,name);

  // Close file
  MPI_File_close(&file);

  // Free buffer
  free(send_buffer);
}
