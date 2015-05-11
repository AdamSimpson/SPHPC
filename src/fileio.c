////////////////////////////////////////////////
// Utility Functions
////////////////////////////////////////////////
#include "fileio.h"

#include <stdio.h>
#include <inttypes.h>

// Write boundary in MPI
void write_MPI(const FluidParticle *const particles,
               const Params *const params,
               const int fileNum)
{
    MPI_File file;
    MPI_Status status;
    char name[64];
    sprintf(name, "/Users/atj/SPHPC/output/sim-%d.bin", fileNum);

    const int num_particles = params->number_fluid_particles_local;

    // How many bytes each process will write
    int rank_write_counts[params->nprocs];
    // alltoall of write counts
    int num_doubles_to_send = 3 * num_particles;
    MPI_Allgather(&num_doubles_to_send, 1, MPI_INT, rank_write_counts, 1, MPI_INT, MPI_COMM_WORLD);
    // Displacement can overflow with int, max size = 8*3*(global num particles)
    uint64_t displacement=0;
    for(int i=0; i<params->rank; i++)
        displacement+=rank_write_counts[i];
    // Displacement is in bytes
    displacement*=sizeof(double);

    // Write x,y,z data to memory
    double *const send_buffer = malloc(num_doubles_to_send*sizeof(double));

    int index=0;
    for(int i=0; i<num_particles; i++) {
        send_buffer[index]   = particles[i].x;
        send_buffer[index+1] = particles[i].y;
        send_buffer[index+2] = particles[i].z;
        index+=3;
    }

    // Open file
    MPI_File_open(MPI_COMM_WORLD, name, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

    // write coordinates to file
    MPI_File_write_at(file, displacement, send_buffer , num_doubles_to_send, MPI_DOUBLE, &status);

    int nbytes;
    MPI_Get_elements(&status, MPI_CHAR, &nbytes);
    printf("rank %d wrote %d bytes(%d particles) to file %s\n",params->rank,nbytes,num_particles,name);

    // Close file
    MPI_File_close(&file);

    // Free buffer
    free(send_buffer);
}
