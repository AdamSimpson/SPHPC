#include "communication.h"
#include "simulation.h"
#include "fluid.h"
#include "debug.h"
#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include "mpi.h"

void AllocateCommunication(struct Communication *const communication) {
  communication->edges.indices_left  = calloc(communication->max_comm_particles, sizeof(int));
  if(communication->edges.indices_left == NULL)
    printf("Could not allocate edge indices left\n");
  communication->edges.indices_right = calloc(communication->max_comm_particles, sizeof(int));
  if(communication->edges.indices_right == NULL)
    printf("Could not allocate edge indices right\n");

  communication->out_of_bounds.indices_left  = calloc(communication->max_comm_particles, sizeof(int));
  if(communication->out_of_bounds.indices_left == NULL)
    printf("Could not allocate oob indices left\n");
  communication->out_of_bounds.indices_right = calloc(communication->max_comm_particles, sizeof(int));
  if(communication->out_of_bounds.indices_right == NULL)
    printf("Could not allocate oob indices right\n");

  // Allocate send and receive buffers
  communication->particle_send_buffer = calloc(communication->max_comm_particles*17, sizeof(double));
  if(communication->particle_send_buffer == NULL)
    printf("Could not allocate send buffer\n");
  communication->particle_recv_buffer = calloc(communication->max_comm_particles*17, sizeof(double));
  if(communication->particle_recv_buffer == NULL)
    printf("Could not allocate recv buffer\n");
}

void FreeCommunication(struct Communication *const communication) {
  free(communication->edges.indices_left);
  free(communication->edges.indices_right);
  free(communication->out_of_bounds.indices_left);
  free(communication->out_of_bounds.indices_right);
  free(communication->particle_send_buffer);
  free(communication->particle_recv_buffer);
}

void InitCommunication(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
}

void FinalizeCommunication() {
  MPI_Finalize();
}

int get_rank() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int get_num_procs() {
  int num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  return num_procs;
}

void PackParticleToBuffer(const struct FluidParticles *const fluid_particles,
                          const int from_index,
                          double *const to_buffer,
                          const int to_index) {
  to_buffer[to_index]      = fluid_particles->x_star[from_index];
  to_buffer[to_index + 1]  = fluid_particles->y_star[from_index];
  to_buffer[to_index + 2]  = fluid_particles->z_star[from_index];
  to_buffer[to_index + 3]  = fluid_particles->x[from_index];
  to_buffer[to_index + 4]  = fluid_particles->y[from_index];
  to_buffer[to_index + 5]  = fluid_particles->z[from_index];
  to_buffer[to_index + 6]  = fluid_particles->v_x[from_index];
  to_buffer[to_index + 7]  = fluid_particles->v_y[from_index];
  to_buffer[to_index + 8]  = fluid_particles->v_z[from_index];
  to_buffer[to_index + 9]  = fluid_particles->dp_x[from_index];
  to_buffer[to_index + 10] = fluid_particles->dp_y[from_index];
  to_buffer[to_index + 11] = fluid_particles->dp_z[from_index];
  to_buffer[to_index + 12] = fluid_particles->w_x[from_index];
  to_buffer[to_index + 13] = fluid_particles->w_y[from_index];
  to_buffer[to_index + 14] = fluid_particles->w_z[from_index];
  to_buffer[to_index + 15] = fluid_particles->density[from_index];
  to_buffer[to_index + 16] = fluid_particles->lambda[from_index];
}

void UnpackBufferToParticle(const double *const from_buffer,
                            const int from_index,
                            struct FluidParticles *const fluid_particles,
                            const int to_index) {
  fluid_particles->x_star[to_index]  = from_buffer[from_index];
  fluid_particles->y_star[to_index]  = from_buffer[from_index + 1];
  fluid_particles->z_star[to_index]  = from_buffer[from_index + 2];
  fluid_particles->x[to_index]       = from_buffer[from_index + 3];
  fluid_particles->y[to_index]       = from_buffer[from_index + 4];
  fluid_particles->z[to_index]       = from_buffer[from_index + 5];
  fluid_particles->v_x[to_index]     = from_buffer[from_index + 6];
  fluid_particles->v_y[to_index]     = from_buffer[from_index + 7];
  fluid_particles->v_z[to_index]     = from_buffer[from_index + 8];
  fluid_particles->dp_x[to_index]    = from_buffer[from_index + 9];
  fluid_particles->dp_y[to_index]    = from_buffer[from_index + 10];
  fluid_particles->dp_z[to_index]    = from_buffer[from_index + 11];
  fluid_particles->w_x[to_index]     = from_buffer[from_index + 12];
  fluid_particles->w_y[to_index]     = from_buffer[from_index + 13];
  fluid_particles->w_z[to_index]     = from_buffer[from_index + 14];
  fluid_particles->density[to_index] = from_buffer[from_index + 15];
  fluid_particles->lambda[to_index]  = from_buffer[from_index + 16];
  fluid_particles->id[to_index] = to_index;
}

// Pack particle struct float components into contiguous memory
void PackHaloComponents(const struct Communication *const communication,
                        const struct FluidParticles *const fluid_particles,
                        double *const packed_send_left,
                        double *const packed_send_right) {

  const struct Edges *const edges = &communication->edges;

  for (int i=0; i<edges->number_particles_left; i++) {
    const int p_index = edges->indices_left[i];
    PackParticleToBuffer(fluid_particles, p_index, packed_send_left, i*17);
  }

  for (int i=0; i<edges->number_particles_right; i++) {
    const int p_index = edges->indices_right[i];
    PackParticleToBuffer(fluid_particles, p_index, packed_send_right, i*17);
  }
}

void UnpackHaloComponents(const struct Params *const params,
                          const double *const packed_recv_left,
                          const double *const packed_recv_right,
                          struct FluidParticles *const fluid_particles) {

  const int num_local = params->number_fluid_particles_local;

  // Unpack halo particles from left rank first
  for (int i=0; i<params->number_halo_particles_left; i++) {
    const int p_index = num_local + i; // "Global" index
    UnpackBufferToParticle(packed_recv_left, i*17, fluid_particles, p_index);
  }

  // Unpack halo particles from right rank second
  for (int i=0; i<params->number_halo_particles_right; i++) {
    const int p_index = num_local
                      + params->number_halo_particles_left + i; // "Global" index
    UnpackBufferToParticle(packed_recv_right, i*17, fluid_particles, p_index);
  }
}

void HaloExchange(struct Communication *const communication,
                  struct Params *const params,
                  struct FluidParticles *const fluid_particles) {

  const int rank = params->rank;
  const int num_procs = params->num_procs;
  struct Edges *const edges = &communication->edges;
  double h = params->smoothing_radius;
  const int num_components = 17;

  // Set edge particle indices and update number
  edges->number_particles_left = 0;
  edges->number_particles_right = 0;
  for (int i=0; i<params->number_fluid_particles_local; ++i) {
    if (fluid_particles->x_star[i] - params->node_start_x <= h)
      edges->indices_left[edges->number_particles_left++] = i;
    else if (params->node_end_x - fluid_particles->x_star[i] <= h)
      edges->indices_right[edges->number_particles_right++] = i;
  }

  int num_moving_left = edges->number_particles_left;
  int num_moving_right = edges->number_particles_right;

  // Setup nodes to left and right of self
  int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  int proc_to_right = (rank == num_procs-1 ? MPI_PROC_NULL : rank+1);

  DEBUG_PRINT("rank %d, halo: will send %d to left, %d to right\n",
              rank, num_moving_left, num_moving_right);

  // Get number of halo particles from right and left
  int num_from_left = 0;
  int num_from_right = 0;
  int tag = 3217;
  // Send number to right and receive from left
  MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag,
               &num_from_left,    1, MPI_INT, proc_to_left,  tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  // Send number to left and receive from right
  tag = 8425;
  MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left, tag,
               &num_from_right,  1, MPI_INT, proc_to_right,tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // Set send/recv buffer points
  double *const packed_send_left  = communication->particle_send_buffer;
  double *const packed_send_right = packed_send_left + num_moving_left*num_components;
  double *const packed_recv_left  = communication->particle_recv_buffer;
  double *const packed_recv_right = packed_recv_left + num_from_left*num_components;

  // Pack halo particle struct components to send
  PackHaloComponents(communication, fluid_particles,
                     packed_send_left, packed_send_right);

  int tagl = 4312;
  int tagr = 5177;

  // Receive halo from left rank
  MPI_Request reqs[4];
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
  MPI_Status statuses[4];
  MPI_Waitall(4, reqs, statuses);

  // Update params struct with halo values
  int num_received_right = 0;
  int num_received_left = 0;
  MPI_Get_count(&statuses[0], MPI_DOUBLE, &num_received_left);
  MPI_Get_count(&statuses[1], MPI_DOUBLE, &num_received_right);
  num_received_left  /= num_components;
  num_received_right /= num_components;

  params->number_halo_particles_left  = num_received_left;
  params->number_halo_particles_right = num_received_right;

  UnpackHaloComponents(params,
                       packed_recv_left,
                       packed_recv_right,
                       fluid_particles);

  DEBUG_PRINT("rank %d, halo: recv %d from left, %d from right\n",
              params->rank,num_received_left,num_received_right);
}

// Pack out of bounds particle components
void PackOOBComponents(const struct Communication *const communication,
                       const struct FluidParticles *const fluid_particles,
                       double *const packed_send_left,
                       double *const packed_send_right) {

  const struct OOB *const oob = &communication->out_of_bounds;

  for (int i=0; i<oob->number_particles_left; i++) {
    const int p_index = oob->indices_left[i];
    PackParticleToBuffer(fluid_particles, p_index, packed_send_left, i*17);
  }

  for (int i=0; i<oob->number_particles_right; i++) {
    const int p_index = oob->indices_right[i];
    PackParticleToBuffer(fluid_particles, p_index, packed_send_right, i*17);
  }
}

void UnpackOOBComponents(const int num_from_left, const int num_from_right,
                         const double *const packed_recv_left,
                         const double *const packed_recv_right,
                         struct Params *const params,
                         struct FluidParticles *const fluid_particles) {

  for (int i=0; i<num_from_left; ++i) {
    const int p_index = params->number_fluid_particles_local;
    UnpackBufferToParticle(packed_recv_left, i*17, fluid_particles, p_index);
    ++params->number_fluid_particles_local;
  }

  for (int i=0; i<num_from_right; ++i) {
    const int p_index = params->number_fluid_particles_local + num_from_left + i;
    UnpackBufferToParticle(packed_recv_left, i*17, fluid_particles, p_index);
    ++params->number_fluid_particles_local;
  }
}

// Transfer particles that are out of node bounds
void TransferOOBParticles(struct Communication *const communication,
                          struct FluidParticles *const particles,
                          struct Params *const params) {

  const int rank = params->rank;
  const int num_procs = params->num_procs;
  struct OOB *const oob = &communication->out_of_bounds;
  const int num_components = 17;

  // Identify out of bound particles
  oob->number_particles_left  = 0;
  oob->number_particles_right = 0;
  for (int i=0; i<params->number_fluid_particles_local; ++i) {
    if (particles->x_star[i] < params->node_start_x && params->rank != 0) {
      oob->indices_left[oob->number_particles_left] = i;
      ++oob->number_particles_left;
    }
    else if (particles->x_star[i] > params->node_end_x &&
       params->rank != params->num_procs-1) {
      oob->indices_right[oob->number_particles_right] = i;
      ++oob->number_particles_right;
    }
  }

  int num_moving_right = oob->number_particles_right;
  int num_moving_left = oob->number_particles_left;

  // Set send buffer points
  double *const packed_send_left  = communication->particle_send_buffer;
  double *const packed_send_right = packed_send_left + (num_moving_left*17);

  PackOOBComponents(communication,
                    particles,
                    packed_send_left,
                    packed_send_right);

  // Move particles from end to fill leaving particles
  for (int i=0; i<oob->number_particles_left; ++i) {
      const int removed_index = oob->indices_left[i];
      MoveParticle(particles, params->number_fluid_particles_local-1, removed_index);
      --params->number_fluid_particles_local;
  }
  for (int i=0; i<oob->number_particles_right; ++i) {
      const int removed_index = oob->indices_right[i];
      MoveParticle(particles, params->number_fluid_particles_local-1, removed_index);
      --params->number_fluid_particles_local;
  }

  // Setup nodes to left and right of self
  int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  int proc_to_right = (rank == num_procs-1 ? MPI_PROC_NULL : rank+1);

  // Get number of particles from right and left
  int num_from_left  = 0;
  int num_from_right = 0;
  int tag = 7006;
  // Send number to right and receive from left
  MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag,
               &num_from_left,    1, MPI_INT,proc_to_left,   tag,
               MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  // Send number to left and receive from right
  tag = 8278;
  MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left,  tag,
               &num_from_right,  1, MPI_INT, proc_to_right, tag,
               MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  // Set recv buffer points
  double *const packed_recv_left  = communication->particle_recv_buffer;
  double *const packed_recv_right = packed_recv_left + (num_from_left*num_components);

  MPI_Status status;
  int num_received_left = 0;
  int num_received_right = 0;

  // Sending to right, recv from left
  tag = 2522;
  MPI_Sendrecv(packed_send_right, num_components*num_moving_right, MPI_DOUBLE, proc_to_right, tag,
               packed_recv_left, num_components*num_from_left, MPI_DOUBLE, proc_to_left, tag,
               MPI_COMM_WORLD, &status);
  MPI_Get_count(&status, MPI_DOUBLE, &num_received_left);
  num_received_left /= num_components;

  // Sending to left, recv from right
  tag = 1165;
  MPI_Sendrecv(packed_send_left, num_components*num_moving_left, MPI_DOUBLE, proc_to_left, tag,
               packed_recv_right, num_components*num_from_right, MPI_DOUBLE, proc_to_right, tag,
               MPI_COMM_WORLD, &status);
  MPI_Get_count(&status, MPI_DOUBLE, &num_received_right);
  num_received_right /= num_components;

  UnpackOOBComponents(num_from_left,num_from_right,
                      packed_recv_left,
                      packed_recv_right,
                      params,
                      particles);

  DEBUG_PRINT("rank %d OOB: sent left %d, right: %d recv left:%d, right: %d\n",
              rank, num_moving_left, num_moving_right, num_received_left, num_received_right);
}

void UpdateHaloLambdas(const struct Communication *const communication,
                       const struct  Params *const params,
                       struct FluidParticles *const fluid_particles) {
  const struct Edges *const edges = &communication->edges;

  const int rank = params->rank;
  const int num_procs = params->num_procs;

  const int num_moving_left  = edges->number_particles_left;
  const int num_moving_right = edges->number_particles_right;

  const int num_from_left  = params->number_halo_particles_left;
  const int num_from_right = params->number_halo_particles_right;

  // Allocate send/recv buffers
  double *const send_lambdas_left  = communication->particle_send_buffer;
  double *const send_lambdas_right = send_lambdas_left + num_moving_left;

  double *const recv_lambdas_left  = communication->particle_recv_buffer;
  double *const recv_lambdas_right = recv_lambdas_left + num_from_left;

  // Pack local halo lambdas
  for (int i=0; i<num_moving_left; i++) {
    const int p_index = edges->indices_left[i];
    send_lambdas_left[i] = fluid_particles->lambda[p_index];
  }
  for (int i=0; i<num_moving_right; i++) {
    const int p_index = edges->indices_right[i];
    send_lambdas_right[i] = fluid_particles->lambda[p_index];
  }

  // Setup nodes to left and right of self
  const int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  const int proc_to_right = (rank == num_procs-1 ? MPI_PROC_NULL : rank+1);

  // Could do async to perhaps increase performance
  // Send lambdas to right and receive from left
  int tag = 784;
  MPI_Sendrecv(send_lambdas_right, num_moving_right, MPI_DOUBLE, proc_to_right, tag,
               recv_lambdas_left, num_from_left, MPI_DOUBLE, proc_to_left, tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // Send lambdas to left and receive from right
  tag = 456;
  MPI_Sendrecv(send_lambdas_left, num_moving_left, MPI_DOUBLE, proc_to_left, tag,
               recv_lambdas_right, num_from_right, MPI_DOUBLE, proc_to_right,tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // Unpack halo particle lambdas
  for (int i=0; i<num_from_left; i++) {
    const int p_index = params->number_fluid_particles_local + i;
    fluid_particles->lambda[p_index] = recv_lambdas_left[i];
  }
  for (int i=0; i<num_from_right; i++) {
    const int p_index = params->number_fluid_particles_local + num_from_left + i;
    fluid_particles->lambda[p_index] = recv_lambdas_right[i];
  }
}

void UpdateHaloPositions(const struct Communication *const communication,
                         const struct Params *const params,
                         struct FluidParticles *const fluid_particles) {
  const struct Edges *const edges = &communication->edges;

  const int rank = params->rank;
  const int num_procs = params->num_procs;

  // x,y,z components required
  const int num_components = 3;
  const int num_moving_left  = num_components*edges->number_particles_left;
  const int num_moving_right = num_components*edges->number_particles_right;

  const int num_from_left  = num_components*params->number_halo_particles_left;
  const int num_from_right = num_components*params->number_halo_particles_right;

  // Set send/recv buffers
  double *const send_positions_left  = communication->particle_send_buffer;
  double *const send_positions_right = send_positions_left + num_moving_left;

  double *const recv_positions_left  = communication->particle_recv_buffer;
  double *const recv_positions_right = recv_positions_left + num_from_left;

  // Pack local edge positions
  for (int i=0; i<num_moving_left; i+=3) {
    const int p_index = edges->indices_left[i/3];
    send_positions_left[i]   = fluid_particles->x_star[p_index];
    send_positions_left[i+1] = fluid_particles->y_star[p_index];
    send_positions_left[i+2] = fluid_particles->z_star[p_index];
  }
  for (int i=0; i<num_moving_right; i+=3) {
    const int p_index = edges->indices_right[i/3];
    send_positions_right[i]   = fluid_particles->x_star[p_index];
    send_positions_right[i+1] = fluid_particles->y_star[p_index];
    send_positions_right[i+2] = fluid_particles->z_star[p_index];
  }

  // Setup nodes to left and right of self
  const int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  const int proc_to_right = (rank == num_procs-1 ? MPI_PROC_NULL : rank+1);

  // Could do async to perhaps increase performance
  // Send positions to right and receive from left
  int tag = 874;
  MPI_Sendrecv(send_positions_right, num_moving_right, MPI_DOUBLE, proc_to_right, tag,
               recv_positions_left,  num_from_left,    MPI_DOUBLE, proc_to_left,  tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // Send positions to left and receive from right
  tag = 546;
  MPI_Sendrecv(send_positions_left,  num_moving_left, MPI_DOUBLE, proc_to_left,  tag,
               recv_positions_right, num_from_right,  MPI_DOUBLE, proc_to_right, tag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Unpack halo particle positions
    for (int i=0; i<num_from_left; i+=3) {
        const int p_index = params->number_fluid_particles_local + i/3;
        fluid_particles->x_star[p_index] = recv_positions_left[i];
        fluid_particles->y_star[p_index] = recv_positions_left[i+1];
        fluid_particles->z_star[p_index] = recv_positions_left[i+2];
    }
    for (int i=0; i<num_from_right; i+=3) {
        const int p_index = params->number_fluid_particles_local
                          + num_from_left/3 + i/3;
        fluid_particles->x_star[p_index] = recv_positions_right[i];
        fluid_particles->y_star[p_index] = recv_positions_right[i+1];
        fluid_particles->z_star[p_index] = recv_positions_right[i+2];
    }
}
