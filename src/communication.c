#include "communication.h"
#include "debug.h"
#include <stdio.h>
#include <stddef.h>
#include <string.h>

void AllocateCommunication(Communication *const communication) {
  // Allocate index arrays
  const size_t index_bytes = communication->max_comm_particles*sizeof(int);

  communication->edges.edge_indices_left  = malloc(index_bytes);
  communication->edges.edge_indices_right = malloc(index_bytes);

  communication->out_of_bounds.oob_indices_left  = malloc(index_bytes);
  communication->out_of_bounds.oob_indices_right = malloc(index_bytes);

  // Allocate send and receive buffers
  const int num_components=3;
  const size_t particle_bytes = communication->max_comm_particles
                               *sizeof(FluidParticle);
  const size_t component_bytes = num_components
                                *communication->max_comm_particles
                                *sizeof(double);

  communication->particle_send_buffer = malloc(particle_bytes);

  communication->halo_components_send_buffer = malloc(component_bytes);
  communication->halo_components_recv_buffer = malloc(component_bytes);
}

void FreeCommunication(Communication *const communication) {
  free(communication->edges.edge_indices_left);
  free(communication->edges.edge_indices_right);
  free(communication->out_of_bounds.oob_indices_left);
  free(communication->out_of_bounds.oob_indices_right);
  free(communication->particle_send_buffer);
  free(communication->halo_components_send_buffer);
  free(communication->halo_components_recv_buffer);
}

void InitCommunication(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  CreateMPITypes();
}

void FinalizeCommunication() {
  FreeMPITypes();
  MPI_Finalize();
}

int get_rank() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

int get_num_procs()
{
  int num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  return num_procs;
}

void CreateMPITypes() {
  MPI_Datatype types[18];
  for (int i=0; i<17; ++i) types[i] = MPI_DOUBLE;
  types[17] = MPI_INT;
  int blocklens[18];
  for (int i=0; i<18; ++i) blocklens[i] = 1;
  MPI_Aint disps[18];
  // Get displacement of each struct member
  disps[0]  = offsetof(FluidParticle, x_star);
  disps[1]  = offsetof(FluidParticle, y_star);
  disps[2]  = offsetof(FluidParticle, z_star);
  disps[3]  = offsetof(FluidParticle, x);
  disps[4]  = offsetof(FluidParticle, y);
  disps[5]  = offsetof(FluidParticle, z);
  disps[6]  = offsetof(FluidParticle, v_x);
  disps[7]  = offsetof(FluidParticle, v_y);
  disps[8]  = offsetof(FluidParticle, v_z);
  disps[9]  = offsetof(FluidParticle, dp_x);
  disps[10] = offsetof(FluidParticle, dp_y);
  disps[11] = offsetof(FluidParticle, dp_z);
  disps[12] = offsetof(FluidParticle, w_x);
  disps[13] = offsetof(FluidParticle, w_y);
  disps[14] = offsetof(FluidParticle, w_z);
  disps[15] = offsetof(FluidParticle, density);
  disps[16] = offsetof(FluidParticle, lambda);
  disps[17] = offsetof(FluidParticle, id);

  MPI_Type_create_struct( 18, blocklens, disps, types, &Particletype );
  MPI_Type_commit( &Particletype );
}

void FreeMPITypes()
{
  MPI_Type_free(&Particletype);
}

void StartHaloExchange(Communication *const communication,
                       const Params *const params,
                       FluidParticle *const fluid_particles) {

  const int rank = params->rank;
  const int num_procs = params->num_procs;
  Edges *const edges = &communication->edges;

  const FluidParticle *p;
  double h = params->smoothing_radius;

  // Set edge particle indices and update number
  edges->number_edge_particles_left = 0;
  edges->number_edge_particles_right = 0;
  for (int i=0; i<params->number_fluid_particles_local; ++i) {
    p = &fluid_particles[i];
    if (p->x_star - params->node_start_x <= h)
      edges->edge_indices_left[edges->number_edge_particles_left++] = i;
    else if (params->node_end_x - p->x_star <= h)
      edges->edge_indices_right[edges->number_edge_particles_right++] = i;
  }

  int num_moving_left = edges->number_edge_particles_left;
  int num_moving_right = edges->number_edge_particles_right;

  // Setup nodes to left and right of self
  const int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  const int proc_to_right = (rank == num_procs-1 ? MPI_PROC_NULL : rank+1);

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

  DEBUG_PRINT("total moving: %d\n", num_moving_right + num_moving_left);

  // Set send buffer points
  FluidParticle *const sendl_buffer = communication->particle_send_buffer;
  FluidParticle *const sendr_buffer = sendl_buffer + num_moving_left;

  // Pack send buffers
  for(int i=0; i<edges->number_edge_particles_left; i++)
    sendl_buffer[i] = fluid_particles[edges->edge_indices_left[i]];
  for(int i=0; i<edges->number_edge_particles_right; i++)
    sendr_buffer[i] = fluid_particles[edges->edge_indices_right[i]];

  //Index to start receiving halo particles
  const int indexToReceiveLeft = params->number_fluid_particles_local;
  const int indexToReceiveRight = indexToReceiveLeft + num_from_left;

  const int tagl = 4312;
  const int tagr = 5177;

  // Receive halo from left rank
  MPI_Irecv(&fluid_particles[indexToReceiveLeft], num_from_left, Particletype,
            proc_to_left, tagl, MPI_COMM_WORLD, &edges->reqs[0]);
  // Receive halo from right rank
  MPI_Irecv(&fluid_particles[indexToReceiveRight], num_from_right, Particletype,
            proc_to_right, tagr, MPI_COMM_WORLD, &edges->reqs[1]);

  // Send halo to right rank
  MPI_Isend(sendr_buffer, num_moving_right, Particletype, proc_to_right, tagl,
            MPI_COMM_WORLD, &edges->reqs[2]);
  // Send halo to left rank
  MPI_Isend(sendl_buffer, num_moving_left, Particletype, proc_to_left, tagr,
            MPI_COMM_WORLD, &edges->reqs[3]);
}

void FinishHaloExchange(Communication *const communication,
                        const FluidParticle *const fluid_particles,
                        Params *const params) {

  // Wait for transfer to complete
  MPI_Status statuses[4];
  MPI_Waitall(4, communication->edges.reqs, statuses);

  int num_received_right = 0;
  int num_received_left = 0;
  MPI_Get_count(&statuses[0], Particletype, &num_received_left);
  MPI_Get_count(&statuses[1], Particletype, &num_received_right);

  params->number_halo_particles_left  = num_received_left;
  params->number_halo_particles_right = num_received_right;

  DEBUG_PRINT("rank %d, halo: recv %d from left, %d from right\n",
              params->rank,num_received_left,num_received_right);
}

// Transfer particles that are out of node bounds
void TransferOOBParticles(const Communication *const communication,
                          FluidParticle *const fluid_particles,
                          Params *const params) {

  const FluidParticle *p;
  const int rank = params->rank;
  const int num_procs = params->num_procs;

  int num_moving_left = 0;
  int num_moving_right = 0;

  // Set send buffers
  FluidParticle *const sendl_buffer = communication->particle_send_buffer;
  FluidParticle *const sendr_buffer = sendl_buffer
                                    + communication->max_comm_particles;

  for (int i=0; i<params->number_fluid_particles_local; i++) {
    p = &fluid_particles[i];

    // Set OOB particle indices and update number
    // remove particle from end of array to fill vacancy
    if (p->x_star < params->node_start_x && params->rank != 0) {
      sendl_buffer[num_moving_left++] = *p;
      fluid_particles[i] = fluid_particles[params->number_fluid_particles_local-1];
      fluid_particles[i].id = i;
      params->number_fluid_particles_local--;
    }
    else if (p->x_star > params->node_end_x &&
             params->rank != params->num_procs-1) {
      sendr_buffer[num_moving_right++] = *p;
      fluid_particles[i] = fluid_particles[params->number_fluid_particles_local-1];
      fluid_particles[i].id = i;
      params->number_fluid_particles_local--;
    }
  }

  // Setup nodes to left and right of self
  const int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  const int proc_to_right = (rank == num_procs-1 ? MPI_PROC_NULL : rank+1);

  // Get number of particles from right and left
  int num_from_left = 0;
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

  // Allocate recieve buffers
  const int num_local = params->number_fluid_particles_local;
  FluidParticle *const recvl_buffer = &fluid_particles[num_local];
  FluidParticle *const recvr_buffer = &fluid_particles[num_local + num_from_left];

  MPI_Status status;

  // Send oob particles to right processor receive oob particles from right processor
  int num_received_left = 0;
  int num_received_right = 0;

  // Sending to right, recv from left
  tag = 2522;
  MPI_Sendrecv(sendr_buffer, num_moving_right, Particletype,proc_to_right, tag,
               recvl_buffer, num_from_left, Particletype, proc_to_left, tag,
               MPI_COMM_WORLD, &status);
  MPI_Get_count(&status, Particletype, &num_received_left);

  // Sending to left, recv from right
  tag = 1165;
  MPI_Sendrecv(sendl_buffer, num_moving_left, Particletype, proc_to_left, tag,
               recvr_buffer, num_from_right, Particletype, proc_to_right, tag,
               MPI_COMM_WORLD, &status);
  MPI_Get_count(&status, Particletype, &num_received_right);

  const int total_received = num_received_right + num_received_left;

  // Update received particle ID
  for (int i=num_local; i<num_local+total_received; i++)
  {
    fluid_particles[i].id = i;
  }

  // Update number of particles
  params->number_fluid_particles_local += total_received;

  DEBUG_PRINT("rank %d OOB: sent left %d, right: %d recv left:%d, right: %d\n",
              rank, num_moving_left, num_moving_right, num_received_left, num_received_right);
}

void UpdateHaloLambdas(const Communication *const communication,
                         const Params *const params,
                         FluidParticle *const fluid_particles)
{
  FluidParticle *p;
  const Edges *const edges = &communication->edges;

  const int rank = params->rank;
  const int num_procs = params->num_procs;

  const int num_moving_left  = edges->number_edge_particles_left;
  const int num_moving_right = edges->number_edge_particles_right;

  const int num_from_left  = params->number_halo_particles_left;
  const int num_from_right = params->number_halo_particles_right;

  // Allocate send/recv buffers
  double *const send_lambdas_left  = communication->halo_components_send_buffer;
  double *const send_lambdas_right = send_lambdas_left + num_moving_left;

  double *const recv_lambdas_left  = communication->halo_components_recv_buffer;
  double *const recv_lambdas_right = recv_lambdas_left + num_from_left;

  // Pack local halo lambdas
  for (int i=0; i<num_moving_left; i++) {
    p = &fluid_particles[edges->edge_indices_left[i]];
    send_lambdas_left[i] = p->lambda;
  }
  for (int i=0; i<num_moving_right; i++) {
    p = &fluid_particles[edges->edge_indices_right[i]];
    send_lambdas_right[i] = p->lambda;
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
    p = &fluid_particles[params->number_fluid_particles_local + i];
    p->lambda = recv_lambdas_left[i];
  }
  for (int i=0; i<num_from_right; i++) {
    p = &fluid_particles[ params->number_fluid_particles_local + num_from_left + i];
    p->lambda = recv_lambdas_right[i];
  }
}

void UpdateHaloPositions(const Communication *const communication,
                         const Params *const params,
                         FluidParticle *const fluid_particles) {
  FluidParticle *p;
  const Edges *const edges = &communication->edges;

  const int rank = params->rank;
  const int num_procs = params->num_procs;

  // x,y,z components required
  const int num_components = 3;
  const int num_moving_left  = num_components*edges->number_edge_particles_left;
  const int num_moving_right = num_components*edges->number_edge_particles_right;

  const int num_from_left  = num_components*params->number_halo_particles_left;
  const int num_from_right = num_components*params->number_halo_particles_right;

  // Allocate send/recv buffers
  // Could combine left/right into single malloc...
  double *const send_positions_left  = communication->halo_components_send_buffer;
  double *const send_positions_right = send_positions_left + num_moving_left;

  double *const recv_positions_left  = communication->halo_components_recv_buffer;
  double *const recv_positions_right = recv_positions_left + num_from_left;

  // Pack local edge positions
  for (int i=0; i<num_moving_left; i+=3) {
    p = &fluid_particles[edges->edge_indices_left[i/3]];
    send_positions_left[i]   = p->x_star;
    send_positions_left[i+1] = p->y_star;
    send_positions_left[i+2] = p->z_star;
  }
  for (int i=0; i<num_moving_right; i+=3) {
    p = &fluid_particles[edges->edge_indices_right[i/3]];
    send_positions_right[i]   = p->x_star;
    send_positions_right[i+1] = p->y_star;
    send_positions_right[i+2] = p->z_star;
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
        p = &fluid_particles[params->number_fluid_particles_local + i/3];;
        p->x_star = recv_positions_left[i];
        p->y_star = recv_positions_left[i+1];
        p->z_star = recv_positions_left[i+2];
    }
    for (int i=0; i<num_from_right; i+=3) {
        p = &fluid_particles[ params->number_fluid_particles_local
                            + num_from_left/3 + i/3];
        p->x_star = recv_positions_right[i];
        p->y_star = recv_positions_right[i+1];
        p->z_star = recv_positions_right[i+2];
    }
}
