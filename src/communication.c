#include "communication.h"
#include "simulation.h"
#include "particles.h"
#include "debug.h"
#include "safe_alloc.h"
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include "mpi.h"

// Number of double components that need to be communicated per particle
static const int num_components = 17;

void InitCommunication(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
}

void AllocateCommunication(struct Communication *const communication) {

  const int num_buffer_indices = communication->max_particles;
  communication->edges.indices_left  = SAFE_ALLOC(num_buffer_indices,
                                                  sizeof(int));
  communication->edges.indices_right = SAFE_ALLOC(num_buffer_indices,
                                                  sizeof(int));
  communication->out_of_bounds.indices_left  = SAFE_ALLOC(num_buffer_indices,
                                                          sizeof(int));
  communication->out_of_bounds.indices_right = SAFE_ALLOC(num_buffer_indices,
                                                          sizeof(int));

  // Allocate send and receive buffers
  const int num_buffer_doubles = communication->max_particles*num_components;
  communication->particle_send_buffer = SAFE_ALLOC(num_buffer_doubles,
                                                   sizeof(double));
  communication->particle_recv_buffer = SAFE_ALLOC(num_buffer_doubles,
                                                   sizeof(double));

}

void FreeCommunication(struct Communication *const communication) {
  free(communication->edges.indices_left);
  free(communication->edges.indices_right);
  free(communication->out_of_bounds.indices_left);
  free(communication->out_of_bounds.indices_right);
  free(communication->particle_send_buffer);
  free(communication->particle_recv_buffer);
}

void FinalizeCommunication() {
  MPI_Finalize();
}

int get_rank() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int get_proc_count() {
  int proc_count;
  MPI_Comm_size(MPI_COMM_WORLD, &proc_count);
  return proc_count;
}

void PackParticleToBuffer(const struct Particles *const particles,
                          const int from_index,
                          double *const to_buffer,
                          const int to_index) {
  to_buffer[to_index]      = particles->x_star[from_index];
  to_buffer[to_index + 1]  = particles->y_star[from_index];
  to_buffer[to_index + 2]  = particles->z_star[from_index];
  to_buffer[to_index + 3]  = particles->x[from_index];
  to_buffer[to_index + 4]  = particles->y[from_index];
  to_buffer[to_index + 5]  = particles->z[from_index];
  to_buffer[to_index + 6]  = particles->v_x[from_index];
  to_buffer[to_index + 7]  = particles->v_y[from_index];
  to_buffer[to_index + 8]  = particles->v_z[from_index];
  to_buffer[to_index + 9]  = particles->dp_x[from_index];
  to_buffer[to_index + 10] = particles->dp_y[from_index];
  to_buffer[to_index + 11] = particles->dp_z[from_index];
  to_buffer[to_index + 12] = particles->w_x[from_index];
  to_buffer[to_index + 13] = particles->w_y[from_index];
  to_buffer[to_index + 14] = particles->w_z[from_index];
  to_buffer[to_index + 15] = particles->density[from_index];
  to_buffer[to_index + 16] = particles->lambda[from_index];
}

void UnpackBufferToParticle(const double *const from_buffer,
                            const int from_index,
                            struct Particles *const particles,
                            const int to_index) {
  particles->x_star[to_index]  = from_buffer[from_index];
  particles->y_star[to_index]  = from_buffer[from_index + 1];
  particles->z_star[to_index]  = from_buffer[from_index + 2];
  particles->x[to_index]       = from_buffer[from_index + 3];
  particles->y[to_index]       = from_buffer[from_index + 4];
  particles->z[to_index]       = from_buffer[from_index + 5];
  particles->v_x[to_index]     = from_buffer[from_index + 6];
  particles->v_y[to_index]     = from_buffer[from_index + 7];
  particles->v_z[to_index]     = from_buffer[from_index + 8];
  particles->dp_x[to_index]    = from_buffer[from_index + 9];
  particles->dp_y[to_index]    = from_buffer[from_index + 10];
  particles->dp_z[to_index]    = from_buffer[from_index + 11];
  particles->w_x[to_index]     = from_buffer[from_index + 12];
  particles->w_y[to_index]     = from_buffer[from_index + 13];
  particles->w_z[to_index]     = from_buffer[from_index + 14];
  particles->density[to_index] = from_buffer[from_index + 15];
  particles->lambda[to_index]  = from_buffer[from_index + 16];
  particles->id[to_index] = to_index;
}

// Pack particle struct float components into contiguous memory
void PackHaloComponents(const struct Communication *const communication,
                        const struct Particles *const particles,
                        double *const packed_send_left,
                        double *const packed_send_right) {

  const struct Edges *const edges = &communication->edges;

  for (int i=0; i<edges->particle_count_left; i++) {
    const int p_index = edges->indices_left[i];
    PackParticleToBuffer(particles, p_index, packed_send_left, i*num_components);
  }

  for (int i=0; i<edges->particle_count_right; i++) {
    const int p_index = edges->indices_right[i];
    PackParticleToBuffer(particles, p_index, packed_send_right, i*num_components);
  }
}

void UnpackHaloComponents(const double *const packed_recv_left,
                          const double *const packed_recv_right,
                          struct Particles *const particles) {

  const int num_local = particles->local_count;

  // Unpack halo particles from left rank first
  for (int i=0; i<particles->halo_count_left; i++) {
    const int p_index = num_local + i; // "Global" index
    UnpackBufferToParticle(packed_recv_left, i*num_components, particles, p_index);
  }

  // Unpack halo particles from right rank second
  for (int i=0; i<particles->halo_count_right; i++) {
    const int p_index = num_local
                      + particles->halo_count_left + i;
    UnpackBufferToParticle(packed_recv_right, i*num_components, particles, p_index);
  }
}

void ExchangeHalo(struct Communication *const communication,
                  struct Params *const params,
                  struct Particles *const particles) {

  const int rank = params->rank;
  const int proc_count = params->proc_count;
  struct Edges *const edges = &communication->edges;
  double h = params->smoothing_radius;

  // Set edge particle indices and update number
  edges->particle_count_left = 0;
  edges->particle_count_right = 0;
  for (int i=0; i<particles->local_count; ++i) {
    if (particles->x_star[i] - params->node_start_x <= h)
      edges->indices_left[edges->particle_count_left++] = i;
    else if (params->node_end_x - particles->x_star[i] <= h)
      edges->indices_right[edges->particle_count_right++] = i;
  }

  int num_moving_left = edges->particle_count_left;
  int num_moving_right = edges->particle_count_right;

  // Setup nodes to left and right of self
  int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  int proc_to_right = (rank == proc_count-1 ? MPI_PROC_NULL : rank+1);

  DEBUG_PRINT("rank %d, halo: will send %d to left, %d to right\n",
              rank, num_moving_left, num_moving_right);

  int tagl = 4312;
  int tagr = 5177;

  // Get number of halo particles from right and left
  int num_from_left = 0;
  int num_from_right = 0;

  // Receive count from left rank
  MPI_Request reqs[4];
  MPI_Irecv(&num_from_left, 1, MPI_INT, proc_to_left,
            tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive count from right rank
  MPI_Irecv(&num_from_right, 1, MPI_INT, proc_to_right,
            tagr, MPI_COMM_WORLD, &reqs[1]);

  // Send count to right rank
  MPI_Isend(&num_moving_right, 1, MPI_INT, proc_to_right,
            tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send count to left rank
  MPI_Isend(&num_moving_left, 1, MPI_INT, proc_to_left,
            tagr, MPI_COMM_WORLD, &reqs[3]);

  // Wait for transfer to complete
  MPI_Status statuses[4];
  MPI_Waitall(4, reqs, statuses);

  if (particles->local_count + num_from_left + num_from_right
      > particles->max_local) {
    printf("Halo overflow: Aborting! \n");
    exit(-1);
  }

  // Set send/recv buffer points
  double *const packed_send_left  = communication->particle_send_buffer;
  double *const packed_send_right = packed_send_left
                                  + num_moving_left*num_components;
  double *const packed_recv_left  = communication->particle_recv_buffer;
  double *const packed_recv_right = packed_recv_left
                                  + num_from_left*num_components;

  // Pack halo particle struct components to send
  PackHaloComponents(communication, particles,
                     packed_send_left, packed_send_right);

  // Receive halo from left rank
  MPI_Irecv(packed_recv_left, num_components*num_from_left, MPI_DOUBLE,
            proc_to_left, tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive halo from right rank
  MPI_Irecv(packed_recv_right, num_components*num_from_right, MPI_DOUBLE,
            proc_to_right, tagr, MPI_COMM_WORLD, &reqs[1]);

  // Send halo to right rank
  MPI_Isend(packed_send_right, num_components*num_moving_right, MPI_DOUBLE,
            proc_to_right, tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send halo to left rank
  MPI_Isend(packed_send_left, num_components*num_moving_left, MPI_DOUBLE,
            proc_to_left, tagr, MPI_COMM_WORLD, &reqs[3]);

  // Wait for transfer to complete
  MPI_Waitall(4, reqs, statuses);

  particles->halo_count_left  = num_from_left;
  particles->halo_count_right = num_from_right;

  UnpackHaloComponents(packed_recv_left,
                       packed_recv_right,
                       particles);

  DEBUG_PRINT("rank %d, halo: recv %d from left, %d from right\n",
              params->rank,num_from_left,num_from_right);
}

// Pack out of bounds particle components
void PackOOBComponents(const struct Communication *const communication,
                       const struct Particles *const particles,
                       double *const packed_send_left,
                       double *const packed_send_right) {

  const struct OOB *const oob = &communication->out_of_bounds;

  for (int i=0; i<oob->particle_count_left; i++) {
    const int p_index = oob->indices_left[i];
    PackParticleToBuffer(particles, p_index, packed_send_left, i*17);
  }

  for (int i=0; i<oob->particle_count_right; i++) {
    const int p_index = oob->indices_right[i];
    PackParticleToBuffer(particles, p_index, packed_send_right, i*17);
  }
}

void UnpackOOBComponents(const int num_from_left, const int num_from_right,
                         const double *const packed_recv_left,
                         const double *const packed_recv_right,
                         struct Particles *const particles) {

  for (int i=0; i<num_from_left; ++i) {
    const int p_index = particles->local_count;
    UnpackBufferToParticle(packed_recv_left, i*17, particles, p_index);
    ++particles->local_count;
  }

  for (int i=0; i<num_from_right; ++i) {
    const int p_index = particles->local_count + num_from_left + i;
    UnpackBufferToParticle(packed_recv_left, i*17, particles, p_index);
    ++particles->local_count;
  }
}

// Transfer particles that are out of node bounds
void ExchangeOOB(struct Communication *const communication,
                          struct Particles *const particles,
                          struct Params *const params) {

  const int rank = params->rank;
  const int proc_count = params->proc_count;
  struct OOB *const oob = &communication->out_of_bounds;

  // Identify out of bound particles
  oob->particle_count_left  = 0;
  oob->particle_count_right = 0;
  for (int i=0; i<particles->local_count; ++i) {
    if (particles->x_star[i] < params->node_start_x && params->rank != 0) {
      oob->indices_left[oob->particle_count_left] = i;
      ++oob->particle_count_left;
    }
    else if (particles->x_star[i] > params->node_end_x &&
       params->rank != params->proc_count-1) {
      oob->indices_right[oob->particle_count_right] = i;
      ++oob->particle_count_right;
    }
  }

  int num_moving_left = oob->particle_count_left;
  int num_moving_right = oob->particle_count_right;

  DEBUG_PRINT("rank %d, OOB: will send %d to left, %d to right\n",
              rank, num_moving_left, num_moving_right);

  // Set send buffer points
  double *const packed_send_left  = communication->particle_send_buffer;
  double *const packed_send_right = packed_send_left
                                  + num_moving_left*num_components;

  PackOOBComponents(communication,
                    particles,
                    packed_send_left,
                    packed_send_right);

  // Move particles from end to fill leaving particles
  for (int i=0; i<oob->particle_count_left; ++i) {
      const int removed_index = oob->indices_left[i];
      CopyParticle(particles, particles->local_count-1, removed_index);
      --particles->local_count;
  }
  for (int i=0; i<oob->particle_count_right; ++i) {
      const int removed_index = oob->indices_right[i];
      CopyParticle(particles, particles->local_count-1, removed_index);
      --particles->local_count;
  }

  // Setup nodes to left and right of self
  int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  int proc_to_right = (rank == proc_count-1 ? MPI_PROC_NULL : rank+1);

  // Get number of particles from right and left
  int num_from_left  = 0;
  int num_from_right = 0;
  int tagl = 7006;
  int tagr = 8278;

  // Receive count from left rank
  MPI_Request reqs[4];
  MPI_Irecv(&num_from_left, 1, MPI_INT, proc_to_left,
            tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive count from right rank
  MPI_Irecv(&num_from_right, 1, MPI_INT, proc_to_right,
            tagr, MPI_COMM_WORLD, &reqs[1]);

  // Send count to right rank
  MPI_Isend(&num_moving_right, 1, MPI_INT, proc_to_right,
            tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send count to left rank
  MPI_Isend(&num_moving_left, 1, MPI_INT, proc_to_left,
            tagr, MPI_COMM_WORLD, &reqs[3]);

  // Wait for transfer to complete
  MPI_Status statuses[4];
  MPI_Waitall(4, reqs, statuses);

  if (particles->local_count + num_from_left + num_from_right
      > particles->max_local) {
    printf("OOB overflow: Aborting! \n");
    exit(-1);
  }

  // Set recv buffer points
  double *const packed_recv_left  = communication->particle_recv_buffer;
  double *const packed_recv_right = packed_recv_left
                                  + num_from_left*num_components;

  // Receive packed doubles from left rank
  MPI_Irecv(packed_recv_left, num_components*num_from_left, MPI_DOUBLE,
            proc_to_left, tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive packed doubles from right rank
  MPI_Irecv(packed_recv_right, num_components*num_from_right, MPI_DOUBLE,
            proc_to_right, tagr, MPI_COMM_WORLD, &reqs[1]);

  // Send packed doubles to right rank
  MPI_Isend(packed_send_right,  num_components*num_moving_right, MPI_DOUBLE,
            proc_to_right, tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send packed doubles to left rank
  MPI_Isend(packed_send_left, num_components*num_moving_left, MPI_DOUBLE,
            proc_to_left, tagr, MPI_COMM_WORLD, &reqs[3]);

  // Wait for transfer to complete
  MPI_Waitall(4, reqs, statuses);

  UnpackOOBComponents(num_from_left, num_from_right,
                      packed_recv_left,
                      packed_recv_right,
                      particles);

  DEBUG_PRINT("rank %d, OOB: recv %d from left, %d from right\n",
              rank, num_from_left, num_from_right);
}

void UpdateHaloLambdas(const struct Communication *const communication,
                       const struct  Params *const params,
                       struct Particles *const particles) {
  const struct Edges *const edges = &communication->edges;

  const int rank = params->rank;
  const int proc_count = params->proc_count;

  const int num_moving_left  = edges->particle_count_left;
  const int num_moving_right = edges->particle_count_right;

  const int num_from_left  = particles->halo_count_left;
  const int num_from_right = particles->halo_count_right;

  // Allocate send/recv buffers
  double *const lambdas_send_left  = communication->particle_send_buffer;
  double *const lambdas_send_right = lambdas_send_left + num_moving_left;

  double *const lambdas_recv_left  = communication->particle_recv_buffer;
  double *const lambdas_recv_right = lambdas_recv_left + num_from_left;

  // Pack local halo lambdas
  for (int i=0; i<num_moving_left; i++) {
    const int p_index = edges->indices_left[i];
    lambdas_send_left[i] = particles->lambda[p_index];
  }
  for (int i=0; i<num_moving_right; i++) {
    const int p_index = edges->indices_right[i];
    lambdas_send_right[i] = particles->lambda[p_index];
  }

  // Setup nodes to left and right of self
  const int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  const int proc_to_right = (rank == proc_count-1 ? MPI_PROC_NULL : rank+1);

  // Send lambdas to right and receive from left
  int tagl = 784;
  int tagr = 456;
  MPI_Request reqs[4];
  // Receive packed doubles from left rank
  MPI_Irecv(lambdas_recv_left, num_from_left, MPI_DOUBLE,
            proc_to_left, tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive packed doubles from right rank
  MPI_Irecv(lambdas_recv_right, num_from_right, MPI_DOUBLE,
            proc_to_right, tagr, MPI_COMM_WORLD, &reqs[1]);

  // Send packed doubles to right rank
  MPI_Isend(lambdas_send_right, num_moving_right, MPI_DOUBLE,
            proc_to_right, tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send packed doubles to left rank
  MPI_Isend(lambdas_send_left, num_moving_left, MPI_DOUBLE,
            proc_to_left, tagr, MPI_COMM_WORLD, &reqs[3]);

  // Wait for transfer to complete
  MPI_Status statuses[4];
  MPI_Waitall(4, reqs, statuses);

  // Unpack halo particle lambdas
  for (int i=0; i<num_from_left; i++) {
    const int p_index = particles->local_count + i;
    particles->lambda[p_index] = lambdas_recv_left[i];
  }
  for (int i=0; i<num_from_right; i++) {
    const int p_index = particles->local_count + num_from_left + i;
    particles->lambda[p_index] = lambdas_recv_right[i];
  }
}

void UpdateHaloPositions(const struct Communication *const communication,
                         const struct Params *const params,
                         struct Particles *const particles) {
  const struct Edges *const edges = &communication->edges;

  const int rank = params->rank;
  const int proc_count = params->proc_count;

  // x,y,z components required
  const int double_components = 3;
  const int num_moving_left  = double_components*edges->particle_count_left;
  const int num_moving_right = double_components*edges->particle_count_right;

  const int num_from_left  = double_components*particles->halo_count_left;
  const int num_from_right = double_components*particles->halo_count_right;

  // Set send/recv buffers
  double *const positions_send_left  = communication->particle_send_buffer;
  double *const positions_send_right = positions_send_left + num_moving_left;

  double *const positions_recv_left  = communication->particle_recv_buffer;
  double *const positions_recv_right = positions_recv_left + num_from_left;

  // Pack local edge positions
  for (int i=0; i<num_moving_left; i+=3) {
    const int p_index = edges->indices_left[i/3];
    positions_send_left[i]   = particles->x_star[p_index];
    positions_send_left[i+1] = particles->y_star[p_index];
    positions_send_left[i+2] = particles->z_star[p_index];
  }
  for (int i=0; i<num_moving_right; i+=3) {
    const int p_index = edges->indices_right[i/3];
    positions_send_right[i]   = particles->x_star[p_index];
    positions_send_right[i+1] = particles->y_star[p_index];
    positions_send_right[i+2] = particles->z_star[p_index];
  }

  // Setup nodes to left and right of self
  const int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  const int proc_to_right = (rank == proc_count-1 ? MPI_PROC_NULL : rank+1);

  // Send lambdas to right and receive from left
  int tagl = 874;
  int tagr = 546;
  MPI_Request reqs[4];
  // Receive packed doubles from left rank
  MPI_Irecv(positions_recv_left, num_from_left, MPI_DOUBLE,
            proc_to_left, tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive packed doubles from right rank
  MPI_Irecv(positions_recv_right, num_from_right, MPI_DOUBLE,
            proc_to_right, tagr, MPI_COMM_WORLD, &reqs[1]);

  // Send packed doubles to right rank
  MPI_Isend(positions_send_right, num_moving_right, MPI_DOUBLE,
            proc_to_right, tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send packed doubles to left rank
  MPI_Isend(positions_send_left, num_moving_left, MPI_DOUBLE,
            proc_to_left, tagr, MPI_COMM_WORLD, &reqs[3]);

  // Wait for transfer to complete
  MPI_Status statuses[4];
  MPI_Waitall(4, reqs, statuses);

  // Unpack halo particle positions
  for (int i=0; i<num_from_left; i+=3) {
    const int p_index = particles->local_count + i/3;
    particles->x_star[p_index] = positions_recv_left[i];
    particles->y_star[p_index] = positions_recv_left[i+1];
    particles->z_star[p_index] = positions_recv_left[i+2];
  }
  for (int i=0; i<num_from_right; i+=3) {
    const int p_index = particles->local_count
                      + num_from_left/3 + i/3;
    particles->x_star[p_index] = positions_recv_right[i];
    particles->y_star[p_index] = positions_recv_right[i+1];
    particles->z_star[p_index] = positions_recv_right[i+2];
  }
}
