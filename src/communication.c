#include "communication.h"
#include "simulation.h"
#include "particles.h"
#include "debug.h"
#include "safe_alloc.h"
#include "thrust_c.h"
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
  communication->send_buffer_left  = SAFE_ALLOC(num_buffer_doubles,
                                                   sizeof(double));
  communication->send_buffer_right = SAFE_ALLOC(num_buffer_doubles,
                                                   sizeof(double));
  communication->recv_buffer_left  = SAFE_ALLOC(num_buffer_doubles,
                                                   sizeof(double));
  communication->recv_buffer_right = SAFE_ALLOC(num_buffer_doubles,
                                                   sizeof(double));
/*
  #pragma acc enter data copyin(                                        \
    communication->edges.indices_left[0:num_buffer_indices],          \
    communication->edges.indices_right[0:num_buffer_indices],         \
    communication->out_of_bounds.indices_left[0:num_buffer_indices],  \
    communication->out_of_bounds.indices_right[0:num_buffer_indices], \
    communication->send_buffer_left[0:num_buffer_doubles],            \
    communication->send_buffer_right[0:num_buffer_doubles],           \
    communication->recv_buffer_left[0:num_buffer_doubles],            \
    communication->recv_buffer_right[0:num_buffer_doubles]            \
  )
*/
}

void FreeCommunication(struct Communication *const communication) {
/*
  #pragma acc exit data delete(                                       \
    communication->edges.indices_left,          \
    communication->edges.indices_right,         \
    communication->out_of_bounds.indices_left,  \
    communication->out_of_bounds.indices_right, \
    communication->send_buffer_left,            \
    communication->send_buffer_right,           \
    communication->recv_buffer_left,            \
    communication->recv_buffer_right            \
  )
*/
  free(communication->edges.indices_left);
  free(communication->edges.indices_right);
  free(communication->out_of_bounds.indices_left);
  free(communication->out_of_bounds.indices_right);
  free(communication->send_buffer_left);
  free(communication->send_buffer_right);
  free(communication->recv_buffer_left);
  free(communication->recv_buffer_right);
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

#pragma acc routine
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

#pragma acc routine
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
                        const struct Particles *const particles) {

  const struct Edges *const edges = &communication->edges;
  double *const packed_send_left  = communication->send_buffer_left;
  double *const packed_send_right = communication->send_buffer_right;

  const int num_particles = particles->local_count;

  #pragma acc parallel loop \
    copyin(particles[:1],                   \
           particles->x_star[:num_particles],\
           particles->y_star[:num_particles],\
           particles->z_star[:num_particles],\
           particles->x[:num_particles],\
           particles->y[:num_particles],\
           particles->z[:num_particles],\
           particles->v_x[:num_particles],\
           particles->v_y[:num_particles],\
           particles->v_z[:num_particles],\
           particles->dp_x[:num_particles],\
           particles->dp_y[:num_particles],\
           particles->dp_z[:num_particles],\
           particles->w_x[:num_particles],\
           particles->w_y[:num_particles],\
           particles->w_z[:num_particles],\
           particles->density[:num_particles],\
           particles->lambda[:num_particles], \
           edges[:1],                          \
           edges->indices_left[:edges->particle_count_left]) \
  copyout(packed_send_left[:edges->particle_count_left*num_components]) default(none)
  for (int i=0; i<edges->particle_count_left; ++i) {
    const int p_index = edges->indices_left[i];
    PackParticleToBuffer(particles, p_index, packed_send_left, i*num_components);
  }

  #pragma acc parallel loop \
    copyin(particles[:1],                   \
           particles->x_star[:num_particles],\
           particles->y_star[:num_particles],\
           particles->z_star[:num_particles],\
           particles->x[:num_particles],\
           particles->y[:num_particles],\
           particles->z[:num_particles],\
           particles->v_x[:num_particles],\
           particles->v_y[:num_particles],\
           particles->v_z[:num_particles],\
           particles->dp_x[:num_particles],\
           particles->dp_y[:num_particles],\
           particles->dp_z[:num_particles],\
           particles->w_x[:num_particles],\
           particles->w_y[:num_particles],\
           particles->w_z[:num_particles],\
           particles->density[:num_particles],\
           particles->lambda[:num_particles], \
           edges[:1],                          \
           edges->indices_right[:edges->particle_count_right]) \
  copyout(packed_send_right[:communication->max_particles*num_components]) default(none)
  for (int i=0; i<edges->particle_count_right; ++i) {
    const int p_index = edges->indices_right[i];
    PackParticleToBuffer(particles, p_index, packed_send_right, i*num_components);
  }
}

void UnpackHaloComponents(const struct Communication *const communication,
                          struct Particles *const particles) {

  const int num_local = particles->local_count;
  const double *const packed_recv_left  = communication->recv_buffer_left;
  const double *const packed_recv_right = communication->recv_buffer_right;

  // Unpack halo particles from left rank first
  int num_particles = particles->local_count + particles->halo_count_left;
  #pragma acc parallel loop \
    copy(particles[:1],                   \
           particles->x_star[:num_particles],\
           particles->y_star[:num_particles],\
           particles->z_star[:num_particles],\
           particles->x[:num_particles],\
           particles->y[:num_particles],\
           particles->z[:num_particles],\
           particles->v_x[:num_particles],\
           particles->v_y[:num_particles],\
           particles->v_z[:num_particles],\
           particles->dp_x[:num_particles],\
           particles->dp_y[:num_particles],\
           particles->dp_z[:num_particles],\
           particles->w_x[:num_particles],\
           particles->w_y[:num_particles],\
           particles->w_z[:num_particles],\
           particles->density[:num_particles],\
           particles->lambda[:num_particles], \
           particles->id[:num_particles]) \
  copyin(packed_recv_left[:particles->halo_count_left*num_components]) default(none)
  for (int i=0; i<particles->halo_count_left; ++i) {
    const int p_index = num_local + i;
    UnpackBufferToParticle(packed_recv_left, i*num_components, particles, p_index);
  }

  // Unpack halo particles from right rank second
  num_particles = particles->local_count + particles->halo_count_left + particles->halo_count_right;
  #pragma acc parallel loop \
    copyin(particles[:1]) \
    copy(particles->x_star[:num_particles],\
           particles->y_star[:num_particles],\
           particles->z_star[:num_particles],\
           particles->x[:num_particles],\
           particles->y[:num_particles],\
           particles->z[:num_particles],\
           particles->v_x[:num_particles],\
           particles->v_y[:num_particles],\
           particles->v_z[:num_particles],\
           particles->dp_x[:num_particles],\
           particles->dp_y[:num_particles],\
           particles->dp_z[:num_particles],\
           particles->w_x[:num_particles],\
           particles->w_y[:num_particles],\
           particles->w_z[:num_particles],\
           particles->density[:num_particles],\
           particles->lambda[:num_particles], \
           particles->id[:num_particles]) \
  copyin(packed_recv_right[:particles->halo_count_right*num_components]) default(none)
  for (int i=0; i<particles->halo_count_right; ++i) {
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

  int num_particles = particles->local_count;
  int max_comm_indices = communication->max_particles;

  #pragma acc data copyin(particles[:1], \
     particles->id[:num_particles],       \
     particles->x_star[:num_particles],   \
     edges[:1])                            \
    copyout(edges->indices_left[:max_comm_indices], \
     edges->indices_right[:max_comm_indices])
  {
  int *ids = particles->id;
  double *x_stars = particles->x_star;
  int *indices_left = edges->indices_left;
  int *indices_right = edges->indices_right;
  #pragma acc host_data use_device(particles)
  {
  // Set edge particle ID's
  CopyIfLessThanOrEqual(h + params->node_start_x,
                        ids,
                        particles->local_count,
                        x_stars,
                        indices_left,
                        &edges->particle_count_left);
  CopyIfGreaterThanOrEqual(params->node_end_x - h,
                           ids,
                           particles->local_count,
                           x_stars,
                           indices_right,
                           &edges->particle_count_right);
  }
  }

  int num_moving_left = edges->particle_count_left;
  int num_moving_right = edges->particle_count_right;

  // Setup nodes to left and right of self
  int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  int proc_to_right = (rank == proc_count-1 ? MPI_PROC_NULL : rank+1);

  DEBUG_PRINT("rank %d, halo: will send %d to %d, %d to %d\n",
              rank, num_moving_left, proc_to_left,
              num_moving_right, proc_to_right);

  // Pack halo particle struct components to send
  PackHaloComponents(communication, particles);

  const int tagl = 4313;
  const int tagr = 5178;
  MPI_Request reqs[4];

  double *const packed_recv_left  = communication->recv_buffer_left;
  double *const packed_recv_right = communication->recv_buffer_right;
  // Receive halo from left rank
  MPI_Irecv(packed_recv_left, num_components*communication->max_particles, MPI_DOUBLE,
            proc_to_left, tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive halo from right rank
  MPI_Irecv(packed_recv_right, num_components*communication->max_particles, MPI_DOUBLE,
            proc_to_right, tagr, MPI_COMM_WORLD, &reqs[1]);

  double *const packed_send_left  = communication->send_buffer_left;
  double *const packed_send_right = communication->send_buffer_right;
  // Send halo to right rank
  MPI_Isend(packed_send_right, num_components*num_moving_right, MPI_DOUBLE,
            proc_to_right, tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send halo to left rank
  MPI_Isend(packed_send_left, num_components*num_moving_left, MPI_DOUBLE,
            proc_to_left, tagr, MPI_COMM_WORLD, &reqs[3]);

  // Wait for transfer to complete
  MPI_Status statuses[4];
  MPI_Waitall(4, reqs, statuses);

  // Need to get number received left and right
  int num_from_left = 0;
  int num_from_right = 0;
  MPI_Get_count(&statuses[0], MPI_DOUBLE, &num_from_left);
  num_from_left /= num_components;
  MPI_Get_count(&statuses[1], MPI_DOUBLE, &num_from_right);
  num_from_right /= num_components;

  particles->halo_count_left  = num_from_left;
  particles->halo_count_right = num_from_right;

  UnpackHaloComponents(communication, particles);

  DEBUG_PRINT("rank %d, halo: recv %d from %d, %d from %d\n",
              params->rank,num_from_left, proc_to_left,
              num_from_right, proc_to_right);
}

// Pack out of bounds particle components
void PackOOBComponents(const struct Communication *const communication,
                       const struct Particles *const particles) {

  const struct OOB *const oob = &communication->out_of_bounds;
  double *const packed_send_left  = communication->send_buffer_left;
  double *const packed_send_right = communication->send_buffer_right;

  int num_particles = particles->local_count;

  #pragma acc parallel loop \
    copyin(particles[:1],                   \
           particles->x_star[:num_particles],\
           particles->y_star[:num_particles],\
           particles->z_star[:num_particles],\
           particles->x[:num_particles],\
           particles->y[:num_particles],\
           particles->z[:num_particles],\
           particles->v_x[:num_particles],\
           particles->v_y[:num_particles],\
           particles->v_z[:num_particles],\
           particles->dp_x[:num_particles],\
           particles->dp_y[:num_particles],\
           particles->dp_z[:num_particles],\
           particles->w_x[:num_particles],\
           particles->w_y[:num_particles],\
           particles->w_z[:num_particles],\
           particles->density[:num_particles],\
           particles->lambda[:num_particles], \
           oob[:1],                          \
           oob->indices_left[:oob->particle_count_left]) \
  copyout(packed_send_left[:oob->particle_count_left*num_components]) default(none)
  for (int i=0; i<oob->particle_count_left; ++i) {
    const int p_index = oob->indices_left[i];
    PackParticleToBuffer(particles, p_index, packed_send_left, i*num_components);
  }

  #pragma acc parallel loop \
    copyin(particles[:1],                   \
           particles->x_star[:num_particles],\
           particles->y_star[:num_particles],\
           particles->z_star[:num_particles],\
           particles->x[:num_particles],\
           particles->y[:num_particles],\
           particles->z[:num_particles],\
           particles->v_x[:num_particles],\
           particles->v_y[:num_particles],\
           particles->v_z[:num_particles],\
           particles->dp_x[:num_particles],\
           particles->dp_y[:num_particles],\
           particles->dp_z[:num_particles],\
           particles->w_x[:num_particles],\
           particles->w_y[:num_particles],\
           particles->w_z[:num_particles],\
           particles->density[:num_particles],\
           particles->lambda[:num_particles], \
           oob[:1],                          \
           oob->indices_right[:oob->particle_count_right]) \
  copyout(packed_send_right[:communication->max_particles*num_components]) default(none)
  for (int i=0; i<oob->particle_count_right; ++i) {
    const int p_index = oob->indices_right[i];
    PackParticleToBuffer(particles, p_index, packed_send_right, i*num_components);
  }
}

void UnpackOOBComponents(const int num_from_left, const int num_from_right,
                         const struct Communication *const communication,
                         struct Particles *const particles) {

  double *const packed_recv_left  = communication->recv_buffer_left;
  double *const packed_recv_right = communication->recv_buffer_right;

  int num_particles = particles->local_count + num_from_left;

  printf("\n right before particles->x_star %p\n", particles->x_star);


  #pragma acc parallel loop \
   copyin(particles[:1])                   \
   copy( particles->x_star[:num_particles],\
          particles->y_star[:num_particles],\
          particles->z_star[:num_particles],\
          particles->x[:num_particles],\
          particles->y[:num_particles],\
          particles->z[:num_particles],\
          particles->v_x[:num_particles],\
          particles->v_y[:num_particles],\
          particles->v_z[:num_particles],\
          particles->dp_x[:num_particles],\
          particles->dp_y[:num_particles],\
          particles->dp_z[:num_particles],\
          particles->w_x[:num_particles],\
          particles->w_y[:num_particles],\
          particles->w_z[:num_particles],\
          particles->density[:num_particles],\
          particles->lambda[:num_particles], \
          particles->id[:num_particles])     \
  copyin(packed_recv_left[:num_from_left*num_components]) default(none)
  for (int i=0; i<num_from_left; ++i) {
    const int p_index = particles->local_count;
    UnpackBufferToParticle(packed_recv_left, i*num_components, particles, p_index);
    ++(particles->local_count);
  }

  num_particles = particles->local_count + num_from_left + num_from_right;

  printf("\n after particles->x_star %p\n", particles->x_star);

  #pragma acc parallel loop \
   copyin(particles[:1])    \
          copy(particles->x_star[:num_particles],\
          particles->y_star[:num_particles],\
          particles->z_star[:num_particles],\
          particles->x[:num_particles],\
          particles->y[:num_particles],\
          particles->z[:num_particles],\
          particles->v_x[:num_particles],\
          particles->v_y[:num_particles],\
          particles->v_z[:num_particles],\
          particles->dp_x[:num_particles],\
          particles->dp_y[:num_particles],\
          particles->dp_z[:num_particles],\
          particles->w_x[:num_particles],\
          particles->w_y[:num_particles],\
          particles->w_z[:num_particles],\
          particles->density[:num_particles],\
          particles->lambda[:num_particles], \
          particles->id[:num_particles])     \
  copyin(packed_recv_right[:num_from_right*num_components]) default(none)
  for (int i=0; i<num_from_right; ++i) {
    const int p_index = particles->local_count;
    UnpackBufferToParticle(packed_recv_right, i*num_components, particles, p_index);
    ++(particles->local_count);
  }

}

// Transfer particles that are out of node bounds
void ExchangeOOB(struct Communication *const communication,
                 struct Particles *const particles,
                 struct Params *const params) {

  const int rank = params->rank;
  const int proc_count = params->proc_count;
  struct OOB *const oob = &communication->out_of_bounds;

  int num_particles = particles->local_count;
  int max_comm_indices = communication->max_particles;

  // use_device doesn't accept struct pointer members so we unpack
  int *ids = particles->id;
  double *x_stars = particles->x_star;
  int *indices_left = oob->indices_left;
  int *indices_right = oob->indices_right;

  #pragma acc data copy(ids[:num_particles]) \
  copyin(x_stars[:num_particles])   \
  copyout(indices_left[:max_comm_indices],  \
    indices_right[:max_comm_indices])
  {

  // Copy OOB particle id's to left node
  #pragma acc host_data use_device(ids, x_stars, indices_left)
  {

  CopyIfLessThan(params->node_start_x,
                 ids,
                 particles->local_count,
                 x_stars,
                 indices_left,
                 &oob->particle_count_left);
  }

  // Copy OOB particle id's to right node
  #pragma acc host_data use_device(ids, x_stars, indices_right)
  {
  CopyIfGreaterThan(params->node_end_x,
                    ids,
                    particles->local_count,
                    x_stars,
                    indices_right,
                    &oob->particle_count_right);
  }

  // Remove OOB particle id's
  #pragma acc host_data use_device(ids, x_stars)
  {
  RemoveIfOutsideBounds(params->node_start_x, params->node_end_x,
                        ids,
                        particles->local_count,
                        x_stars);
  }

  }

  // Pack removed particle components
  PackOOBComponents(communication,
                    particles);

  // Set new particle count
  int num_moving_left = oob->particle_count_left;
  int num_moving_right = oob->particle_count_right;
  particles->local_count -= (num_moving_left + num_moving_right);

  // Use ID component to reorganize particles array
  for(int i=0; i<particles->local_count; ++i) {
      const int current_id = particles->id[i];
      CopyParticle(particles, current_id, i);
  }

  // Setup nodes to left and right of self
  int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  int proc_to_right = (rank == proc_count-1 ? MPI_PROC_NULL : rank+1);

  DEBUG_PRINT("rank %d, OOB: will send %d to %d, %d to %d\n",
              rank, num_moving_left, proc_to_left,
              num_moving_right, proc_to_right);

  int tagl = 7007;
  int tagr = 8279;
  MPI_Request reqs[4];

  double *const packed_recv_left  = communication->recv_buffer_left;
  double *const packed_recv_right = communication->recv_buffer_right;
  // Receive packed doubles from left rank
  MPI_Irecv(packed_recv_left, num_components*communication->max_particles, MPI_DOUBLE,
            proc_to_left, tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive packed doubles from right rank
  MPI_Irecv(packed_recv_right, num_components*communication->max_particles, MPI_DOUBLE,
            proc_to_right, tagr, MPI_COMM_WORLD, &reqs[1]);

  double *const packed_send_left  = communication->send_buffer_left;
  double *const packed_send_right = communication->send_buffer_right;
  // Send packed doubles to right rank
  MPI_Isend(packed_send_right,  num_components*num_moving_right, MPI_DOUBLE,
            proc_to_right, tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send packed doubles to left rank
  MPI_Isend(packed_send_left, num_components*num_moving_left, MPI_DOUBLE,
            proc_to_left, tagr, MPI_COMM_WORLD, &reqs[3]);

  // Wait for transfer to complete
  MPI_Status statuses[4];
  MPI_Waitall(4, reqs, statuses);

  // Need to get number received left and right
  int num_from_left = 0;
  int num_from_right = 0;
  MPI_Get_count(&statuses[0], MPI_DOUBLE, &num_from_left);
  num_from_left /= num_components;
  MPI_Get_count(&statuses[1], MPI_DOUBLE, &num_from_right);
  num_from_right /= num_components;

  printf("\nbefore particles->x_star %p\n", particles->x_star);

  UnpackOOBComponents(num_from_left, num_from_right,
                      communication,
                      particles);

  DEBUG_PRINT("rank %d, OOB: recv %d from %d, %d from %d: count :%d\n",
              rank, num_from_left, proc_to_left, num_from_right, proc_to_right,
              particles->local_count);
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

  double *const lambdas_send_left  = communication->send_buffer_left;
  double *const lambdas_send_right = communication->send_buffer_right;

  // Pack local halo lambdas
  for (int i=0; i<num_moving_left; ++i) {
    const int p_index = edges->indices_left[i];
    lambdas_send_left[i] = particles->lambda[p_index];
  }
  for (int i=0; i<num_moving_right; ++i) {
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

  double *const lambdas_recv_left  = communication->recv_buffer_left;
  double *const lambdas_recv_right = communication->recv_buffer_right;
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
  for (int i=0; i<num_from_left; ++i) {
    const int p_index = particles->local_count + i;
    particles->lambda[p_index] = lambdas_recv_left[i];
  }
  for (int i=0; i<num_from_right; ++i) {
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
  double *const positions_send_left  = communication->send_buffer_left;
  double *const positions_send_right = communication->send_buffer_right;

  double *const positions_recv_left  = communication->recv_buffer_left;
  double *const positions_recv_right = communication->recv_buffer_right;

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
