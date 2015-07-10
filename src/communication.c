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

static const int g_num_components = 9;

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
  const int num_buffer_doubles = communication->max_particles*g_num_components;
  communication->send_buffer_left  = SAFE_ALLOC(num_buffer_doubles,
                                                   sizeof(double));
  communication->send_buffer_right = SAFE_ALLOC(num_buffer_doubles,
                                                   sizeof(double));
  communication->recv_buffer_left  = SAFE_ALLOC(num_buffer_doubles,
                                                   sizeof(double));
  communication->recv_buffer_right = SAFE_ALLOC(num_buffer_doubles,
                                                   sizeof(double));

  #pragma acc enter data copyin( communication[:1],                   \
    communication->edges.indices_left[0:num_buffer_indices],          \
    communication->edges.indices_right[0:num_buffer_indices],         \
    communication->out_of_bounds.indices_left[0:num_buffer_indices],  \
    communication->out_of_bounds.indices_right[0:num_buffer_indices], \
    communication->send_buffer_left[0:num_buffer_doubles],            \
    communication->send_buffer_right[0:num_buffer_doubles],           \
    communication->recv_buffer_left[0:num_buffer_doubles],            \
    communication->recv_buffer_right[0:num_buffer_doubles]            \
  )
}

void FreeCommunication(struct Communication *const communication) {

  #pragma acc exit data delete(                 \
    communication->edges.indices_left,          \
    communication->edges.indices_right,         \
    communication->out_of_bounds.indices_left,  \
    communication->out_of_bounds.indices_right, \
    communication->send_buffer_left,            \
    communication->send_buffer_right,           \
    communication->recv_buffer_left,            \
    communication->recv_buffer_right            \
  )

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

// Pack particle struct float components into contiguous memory
void PackHaloComponents(const struct Communication *const communication,
                        const struct Particles *const particles) {

  const struct Edges *const edges = &communication->edges;
  double *const packed_send_left  = communication->send_buffer_left;
  double *const packed_send_right = communication->send_buffer_right;

  const int num_particles = particles->local_count;
  const int count_left = edges->particle_count_left;
  const int count_right = edges->particle_count_right;

  // OpenACC can't reliably handle SoA...
  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;
  double *x      = particles->x;
  double *y      = particles->y;
  double *z      = particles->z;
  double *v_x    = particles->v_x;
  double *v_y    = particles->v_y;
  double *v_z    = particles->v_z;

  #pragma acc parallel loop present(particles, \
                                    x_star, y_star, z_star, \
                                    x, y, z, \
                                    v_x, v_y, v_z, \
                                    edges, edges->indices_left, packed_send_left) default(none)
  for (int i=0; i<count_left; ++i) {
    const int p_index = edges->indices_left[i];
    packed_send_left[i*g_num_components]      = x_star[p_index];
    packed_send_left[i*g_num_components + 1]  = y_star[p_index];
    packed_send_left[i*g_num_components + 2]  = z_star[p_index];
    packed_send_left[i*g_num_components + 3]  = x[p_index];
    packed_send_left[i*g_num_components + 4]  = y[p_index];
    packed_send_left[i*g_num_components + 5]  = z[p_index];
    packed_send_left[i*g_num_components + 6]  = v_x[p_index];
    packed_send_left[i*g_num_components + 7]  = v_y[p_index];
    packed_send_left[i*g_num_components + 8]  = v_z[p_index];
  }

  #pragma acc parallel loop present(particles, \
                                    x_star, y_star, z_star, \
                                    x, y, z, \
                                    v_x, v_y, v_z, \
                                    edges, edges->indices_right, packed_send_right) default(none)
  for (int i=0; i<count_right; ++i) {
    const int p_index = edges->indices_right[i];
    packed_send_right[i*g_num_components]      = x_star[p_index];
    packed_send_right[i*g_num_components + 1]  = y_star[p_index];
    packed_send_right[i*g_num_components + 2]  = z_star[p_index];
    packed_send_right[i*g_num_components + 3]  = x[p_index];
    packed_send_right[i*g_num_components + 4]  = y[p_index];
    packed_send_right[i*g_num_components + 5]  = z[p_index];
    packed_send_right[i*g_num_components + 6]  = v_x[p_index];
    packed_send_right[i*g_num_components + 7]  = v_y[p_index];
    packed_send_right[i*g_num_components + 8]  = v_z[p_index];
  }
}

void UnpackHaloComponents(const struct Communication *const communication,
                          struct Particles *const particles) {

  const int num_local = particles->local_count;
  const double *const packed_recv_left  = communication->recv_buffer_left;
  const double *const packed_recv_right = communication->recv_buffer_right;

  const int count_left = particles->halo_count_left;
  const int count_right = particles->halo_count_right;

  // OpenACC can't reliably handle SoA...
  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;
  double *x      = particles->x;
  double *y      = particles->y;
  double *z      = particles->z;
  double *v_x    = particles->v_x;
  double *v_y    = particles->v_y;
  double *v_z    = particles->v_z;
  int *id         = particles->id;

  // Unpack halo particles from left rank first
  #pragma acc parallel loop present(particles, \
                                    x_star, y_star, z_star, \
                                    x, y, z, \
                                    v_x, v_y, v_z, \
                                    id, \
                                    packed_recv_left) default(none)
  for (int i=0; i<count_left; ++i) {
    const int p_index = num_local + i;
    x_star[p_index]  = packed_recv_left[i*g_num_components];
    y_star[p_index]  = packed_recv_left[i*g_num_components + 1];
    z_star[p_index]  = packed_recv_left[i*g_num_components + 2];
    x[p_index]       = packed_recv_left[i*g_num_components + 3];
    y[p_index]       = packed_recv_left[i*g_num_components + 4];
    z[p_index]       = packed_recv_left[i*g_num_components + 5];
    v_x[p_index]     = packed_recv_left[i*g_num_components + 6];
    v_y[p_index]     = packed_recv_left[i*g_num_components + 7];
    v_z[p_index]     = packed_recv_left[i*g_num_components + 8];
    id[p_index] = p_index;
  }

  // Unpack halo particles from right rank second
  #pragma acc parallel loop present(particles, \
                                    x_star, y_star, z_star, \
                                    x, y, z, \
                                    v_x, v_y, v_z, \
                                    id, \
                                    packed_recv_right) default(none)
  for (int i=0; i<count_right; ++i) {
    const int p_index = num_local + count_left + i;
    x_star[p_index]  = packed_recv_right[i*g_num_components];
    y_star[p_index]  = packed_recv_right[i*g_num_components + 1];
    z_star[p_index]  = packed_recv_right[i*g_num_components + 2];
    x[p_index]       = packed_recv_right[i*g_num_components + 3];
    y[p_index]       = packed_recv_right[i*g_num_components + 4];
    z[p_index]       = packed_recv_right[i*g_num_components + 5];
    v_x[p_index]     = packed_recv_right[i*g_num_components + 6];
    v_y[p_index]     = packed_recv_right[i*g_num_components + 7];
    v_z[p_index]     = packed_recv_right[i*g_num_components + 8];
    id[p_index] = p_index;
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

  int *ids = particles->id;
  double *x_stars = particles->x_star;
  int *indices_left = edges->indices_left;
  int *indices_right = edges->indices_right;

  #pragma acc host_data use_device(ids, x_stars, indices_left, indices_right)
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

  #pragma acc host_data use_device(packed_recv_left, packed_recv_right)
  {
  // Receive halo from left rank
  MPI_Irecv(packed_recv_left, g_num_components*communication->max_particles, MPI_DOUBLE,
            proc_to_left, tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive halo from right rank
  MPI_Irecv(packed_recv_right, g_num_components*communication->max_particles, MPI_DOUBLE,
            proc_to_right, tagr, MPI_COMM_WORLD, &reqs[1]);
  }

  double *const packed_send_left  = communication->send_buffer_left;
  double *const packed_send_right = communication->send_buffer_right;

  #pragma acc host_data use_device(packed_send_left, packed_send_right)
  {
  // Send halo to right rank
  MPI_Isend(packed_send_right, g_num_components*num_moving_right, MPI_DOUBLE,
            proc_to_right, tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send halo to left rank
  MPI_Isend(packed_send_left, g_num_components*num_moving_left, MPI_DOUBLE,
            proc_to_left, tagr, MPI_COMM_WORLD, &reqs[3]);
  }

  // Wait for transfer to complete
  MPI_Status statuses[4];
  MPI_Waitall(4, reqs, statuses);

  // Need to get number received left and right
  int num_from_left = 0;
  int num_from_right = 0;
  MPI_Get_count(&statuses[0], MPI_DOUBLE, &num_from_left);
  num_from_left /= g_num_components;
  MPI_Get_count(&statuses[1], MPI_DOUBLE, &num_from_right);
  num_from_right /= g_num_components;

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
  int count_left = oob->particle_count_left;
  int count_right = oob->particle_count_right;
  int comp_count_left = count_left*g_num_components;
  int comp_count_right = count_right*g_num_components;

  // OpenACC can't reliably handle SoA...
  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;
  double *x      = particles->x;
  double *y      = particles->y;
  double *z      = particles->z;
  double *v_x    = particles->v_x;
  double *v_y    = particles->v_y;
  double *v_z    = particles->v_z;

  #pragma acc parallel loop present(particles, \
                                      x_star, y_star, z_star, \
                                      x, y, z, \
                                      v_x, v_y, v_z, \
                                    oob, oob->indices_left, packed_send_left) default(none)
  for (int i=0; i<count_left; ++i) {
    const int p_index = oob->indices_left[i];
    packed_send_left[i*g_num_components]      = x_star[p_index];
    packed_send_left[i*g_num_components + 1]  = y_star[p_index];
    packed_send_left[i*g_num_components + 2]  = z_star[p_index];
    packed_send_left[i*g_num_components + 3]  = x[p_index];
    packed_send_left[i*g_num_components + 4]  = y[p_index];
    packed_send_left[i*g_num_components + 5]  = z[p_index];
    packed_send_left[i*g_num_components + 6]  = v_x[p_index];
    packed_send_left[i*g_num_components + 7]  = v_y[p_index];
    packed_send_left[i*g_num_components + 8]  = v_z[p_index];
  }

  #pragma acc parallel loop present(particles, \
                                      x_star, y_star, z_star, \
                                      x, y, z, \
                                      v_x, v_y, v_z, \
                                    oob, oob->indices_right, packed_send_right) default(none)
  for (int i=0; i<count_right; ++i) {
    const int p_index = oob->indices_right[i];
    packed_send_right[i*g_num_components]      = x_star[p_index];
    packed_send_right[i*g_num_components + 1]  = y_star[p_index];
    packed_send_right[i*g_num_components + 2]  = z_star[p_index];
    packed_send_right[i*g_num_components + 3]  = x[p_index];
    packed_send_right[i*g_num_components + 4]  = y[p_index];
    packed_send_right[i*g_num_components + 5]  = z[p_index];
    packed_send_right[i*g_num_components + 6]  = v_x[p_index];
    packed_send_right[i*g_num_components + 7]  = v_y[p_index];
    packed_send_right[i*g_num_components + 8]  = v_z[p_index];
  }
}

// Tomfoolery because particles->local_count is kept updated on device.
void UnpackOOBComponents(const int num_from_left, const int num_from_right,
                         const struct Communication *const communication,
                         struct Particles *const particles) {

  double *const packed_recv_left  = communication->recv_buffer_left;
  double *const packed_recv_right = communication->recv_buffer_right;

  int local_count = particles->local_count;

  // OpenACC can't reliably handle SoA...
  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;
  double *x      = particles->x;
  double *y      = particles->y;
  double *z      = particles->z;
  double *v_x    = particles->v_x;
  double *v_y    = particles->v_y;
  double *v_z    = particles->v_z;
  int *id         = particles->id;

  #pragma acc parallel loop present(particles, \
                                    x_star, y_star, z_star, \
                                    x, y, z, \
                                    v_x, v_y, v_z, \
                                    id, \
                                    packed_recv_left) default(none)
  for (int i=0; i<num_from_left; ++i) {
    const int p_index = local_count + i;
    x_star[p_index]  = packed_recv_left[i*g_num_components];
    y_star[p_index]  = packed_recv_left[i*g_num_components + 1];
    z_star[p_index]  = packed_recv_left[i*g_num_components + 2];
    x[p_index]       = packed_recv_left[i*g_num_components + 3];
    y[p_index]       = packed_recv_left[i*g_num_components + 4];
    z[p_index]       = packed_recv_left[i*g_num_components + 5];
    v_x[p_index]     = packed_recv_left[i*g_num_components + 6];
    v_y[p_index]     = packed_recv_left[i*g_num_components + 7];
    v_z[p_index]     = packed_recv_left[i*g_num_components + 8];
    id[p_index] = p_index;
  }

  #pragma acc parallel loop present(particles, \
                                    x_star, y_star, z_star, \
                                    x, y, z, \
                                    v_x, v_y, v_z, \
                                    id, \
                                    packed_recv_right) default(none)
  for (int i=0; i<num_from_right; ++i) {
    const int p_index = local_count + num_from_left + i;
    x_star[p_index]  = packed_recv_right[i*g_num_components];
    y_star[p_index]  = packed_recv_right[i*g_num_components + 1];
    z_star[p_index]  = packed_recv_right[i*g_num_components + 2];
    x[p_index]       = packed_recv_right[i*g_num_components + 3];
    y[p_index]       = packed_recv_right[i*g_num_components + 4];
    z[p_index]       = packed_recv_right[i*g_num_components + 5];
    v_x[p_index]     = packed_recv_right[i*g_num_components + 6];
    v_y[p_index]     = packed_recv_right[i*g_num_components + 7];
    v_z[p_index]     = packed_recv_right[i*g_num_components + 8];
    id[p_index] = p_index;
  }

  particles->local_count += (num_from_left + num_from_right);
}

// Transfer particles that are out of node bounds
void ExchangeOOB(struct Communication *const communication,
                 struct Particles *const particles,
                 struct Params *const params) {

  const int rank = params->rank;
  const int proc_count = params->proc_count;
  struct OOB *const oob = &communication->out_of_bounds;

  // OpenACC can't reliably handle SoA...
  double *x_star = particles->x_star;
  double *y_star = particles->y_star;
  double *z_star = particles->z_star;
  double *x    = particles->x;
  double *y    = particles->y;
  double *z    = particles->z;
  double *v_x    = particles->v_x;
  double *v_y    = particles->v_y;
  double *v_z    = particles->v_z;
  int *id         = particles->id;
  int *indices_left = oob->indices_left;
  int *indices_right = oob->indices_right;

  #pragma acc host_data use_device(id, x_star, indices_left, indices_right)
  {
  // Copy OOB particle id's to left node
  CopyIfLessThan(params->node_start_x,
                 id,
                 particles->local_count,
                 x_star,
                 indices_left,
                 &oob->particle_count_left);

  // Copy OOB particle id's to right node
  CopyIfGreaterThan(params->node_end_x,
                    id,
                    particles->local_count,
                    x_star,
                    indices_right,
                    &oob->particle_count_right);
  }

  // Pack removed particle components
  PackOOBComponents(communication, particles);

  #pragma acc host_data use_device(id, x_star, y_star, z_star,  \
                                   x, y, z, v_x, v_y, v_z)
  {
  // Particles->id not correct at this point, but not used...should just remove particles->id
  RemoveIfOutsideBounds3(params->node_start_x, params->node_end_x,
                         particles->local_count,
                         x_star, y_star, z_star,
                         x, y, z,
                         v_x, v_y, v_z);

  }

  // Set new particle count
  int num_moving_left = oob->particle_count_left;
  int num_moving_right = oob->particle_count_right;

  particles->local_count -= (num_moving_left + num_moving_right);

  // Must reset id's now
  // Should get rid of particle ID and use counting iterator(?)
  const int local_count = particles->local_count;
  #pragma acc parallel loop present(id)
  for(int i=0; i<local_count; i++) {
    id[i] = i;
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

  #pragma acc host_data use_device(packed_recv_left, packed_recv_right)
  {
  // Receive packed doubles from left rank
  MPI_Irecv(packed_recv_left, g_num_components*communication->max_particles, MPI_DOUBLE,
            proc_to_left, tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive packed doubles from right rank
  MPI_Irecv(packed_recv_right, g_num_components*communication->max_particles, MPI_DOUBLE,
            proc_to_right, tagr, MPI_COMM_WORLD, &reqs[1]);
  }

  double *const packed_send_left  = communication->send_buffer_left;
  double *const packed_send_right = communication->send_buffer_right;

  #pragma acc host_data use_device(packed_send_left, packed_send_right)
  {
  // Send packed doubles to right rank
  MPI_Isend(packed_send_right,  g_num_components*num_moving_right, MPI_DOUBLE,
            proc_to_right, tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send packed doubles to left rank
  MPI_Isend(packed_send_left, g_num_components*num_moving_left, MPI_DOUBLE,
            proc_to_left, tagr, MPI_COMM_WORLD, &reqs[3]);
  }

  // Wait for transfer to complete
  MPI_Status statuses[4];
  MPI_Waitall(4, reqs, statuses);

  // Need to get number received left and right
  int num_from_left = 0;
  int num_from_right = 0;
  MPI_Get_count(&statuses[0], MPI_DOUBLE, &num_from_left);
  num_from_left /= g_num_components;
  MPI_Get_count(&statuses[1], MPI_DOUBLE, &num_from_right);
  num_from_right /= g_num_components;

  UnpackOOBComponents(num_from_left, num_from_right,
                      communication,
                      particles);

  DEBUG_PRINT("rank %d, OOB: recv %d from %d, %d from %d: count :%d\n",
              rank, num_from_left, proc_to_left, num_from_right, proc_to_right,
              particles->local_count);

}

void UpdateHaloScalar(const struct Communication *const communication,
                      const struct  Params *const params,
                      struct Particles *const particles,
                      double *scalars) {
  const struct Edges *const edges = &communication->edges;

  const int rank = params->rank;
  const int proc_count = params->proc_count;

  const int num_moving_left  = edges->particle_count_left;
  const int num_moving_right = edges->particle_count_right;

  const int num_from_left  = particles->halo_count_left;
  const int num_from_right = particles->halo_count_right;

  double *const scalars_send_left  = communication->send_buffer_left;
  double *const scalars_send_right = communication->send_buffer_right;

  // Pack local halo scalars
  #pragma acc parallel loop present(edges, edges->indices_left, scalars_send_left, scalars)
  for (int i=0; i<num_moving_left; ++i) {
    const int p_index = edges->indices_left[i];
    scalars_send_left[i] = scalars[p_index];
  }
  #pragma acc parallel loop present(edges, edges->indices_right, scalars_send_right, scalars)
  for (int i=0; i<num_moving_right; ++i) {
    const int p_index = edges->indices_right[i];
    scalars_send_right[i] = scalars[p_index];
  }

  // Setup nodes to left and right of self
  const int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  const int proc_to_right = (rank == proc_count-1 ? MPI_PROC_NULL : rank+1);

  // Send scalars to right and receive from left
  int tagl = 784;
  int tagr = 456;
  MPI_Request reqs[4];

  double *const scalars_recv_left  = communication->recv_buffer_left;
  double *const scalars_recv_right = communication->recv_buffer_right;

  #pragma acc host_data use_device(scalars_recv_left, scalars_recv_right, scalars_send_left, scalars_send_right)
  {
  // Receive packed doubles from left rank
  MPI_Irecv(scalars_recv_left, num_from_left, MPI_DOUBLE,
            proc_to_left, tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive packed doubles from right rank
  MPI_Irecv(scalars_recv_right, num_from_right, MPI_DOUBLE,
            proc_to_right, tagr, MPI_COMM_WORLD, &reqs[1]);

  // Send packed doubles to right rank
  MPI_Isend(scalars_send_right, num_moving_right, MPI_DOUBLE,
            proc_to_right, tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send packed doubles to left rank
  MPI_Isend(scalars_send_left, num_moving_left, MPI_DOUBLE,
            proc_to_left, tagr, MPI_COMM_WORLD, &reqs[3]);
  }

  // Wait for transfer to complete
  MPI_Status statuses[4];
  MPI_Waitall(4, reqs, statuses);

  const int local_count = particles->local_count;

  // Unpack halo particle scalars
  #pragma acc parallel loop present(scalars_recv_left, scalars)
  for (int i=0; i<num_from_left; ++i) {
    const int p_index = local_count + i;
    scalars[p_index] = scalars_recv_left[i];
  }
  #pragma acc parallel loop present(scalars_recv_right, scalars)
  for (int i=0; i<num_from_right; ++i) {
    const int p_index = local_count + num_from_left + i;
    scalars[p_index] = scalars_recv_right[i];
  }
}

void UpdateHaloTuple(const struct Communication *const communication,
                     const struct Params *const params,
                     struct Particles *const particles,
                     double *const t1, double *const t2, double *const t3) {
  const struct Edges *const edges = &communication->edges;

  const int rank = params->rank;
  const int proc_count = params->proc_count;

  // tuple(t1, t2, t3) components required
  const int double_components = 3;
  const int num_moving_left  = double_components*edges->particle_count_left;
  const int num_moving_right = double_components*edges->particle_count_right;

  const int num_from_left  = double_components*particles->halo_count_left;
  const int num_from_right = double_components*particles->halo_count_right;

  // Set send/recv buffers
  double *const tuples_send_left  = communication->send_buffer_left;
  double *const tuples_send_right = communication->send_buffer_right;

  double *const tuples_recv_left  = communication->recv_buffer_left;
  double *const tuples_recv_right = communication->recv_buffer_right;

  // Pack local edge tuples
  #pragma acc parallel loop present(edges, edges->indices_left, tuples_send_left, t1, t2, t3)
  for (int i=0; i<num_moving_left; i+=3) {
    const int p_index = edges->indices_left[i/3];
    tuples_send_left[i]   = t1[p_index];
    tuples_send_left[i+1] = t2[p_index];
    tuples_send_left[i+2] = t3[p_index];
  }
  #pragma acc parallel loop present(edges, edges->indices_right, tuples_send_right, t1, t2, t3)
  for (int i=0; i<num_moving_right; i+=3) {
    const int p_index = edges->indices_right[i/3];
    tuples_send_right[i]   = t1[p_index];
    tuples_send_right[i+1] = t2[p_index];
    tuples_send_right[i+2] = t3[p_index];
  }

  // Setup nodes to left and right of self
  const int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
  const int proc_to_right = (rank == proc_count-1 ? MPI_PROC_NULL : rank+1);

  // Send tuples to right and receive from left
  int tagl = 874;
  int tagr = 546;
  MPI_Request reqs[4];

  #pragma acc host_data use_device(tuples_recv_left, tuples_recv_right, tuples_send_left, tuples_send_right)
  {
  // Receive packed doubles from left rank
  MPI_Irecv(tuples_recv_left, num_from_left, MPI_DOUBLE,
            proc_to_left, tagl, MPI_COMM_WORLD, &reqs[0]);
  // Receive packed doubles from right rank
  MPI_Irecv(tuples_recv_right, num_from_right, MPI_DOUBLE,
            proc_to_right, tagr, MPI_COMM_WORLD, &reqs[1]);

  // Send packed doubles to right rank
  MPI_Isend(tuples_send_right, num_moving_right, MPI_DOUBLE,
            proc_to_right, tagl, MPI_COMM_WORLD, &reqs[2]);
  // Send packed doubles to left rank
  MPI_Isend(tuples_send_left, num_moving_left, MPI_DOUBLE,
            proc_to_left, tagr, MPI_COMM_WORLD, &reqs[3]);
  }

  // Wait for transfer to complete
  MPI_Status statuses[4];
  MPI_Waitall(4, reqs, statuses);

  const int local_count = particles->local_count;

  // Unpack halo particle tuples
  #pragma acc parallel loop present(tuples_recv_left, t1, t2, t3)
  for (int i=0; i<num_from_left; i+=3) {
    const int p_index = local_count + i/3;
    t1[p_index] = tuples_recv_left[i];
    t2[p_index] = tuples_recv_left[i+1];
    t3[p_index] = tuples_recv_left[i+2];
  }
  #pragma acc parallel loop present(tuples_recv_right, t1, t2, t3)
  for (int i=0; i<num_from_right; i+=3) {
    const int p_index = local_count
                      + num_from_left/3 + i/3;
    t1[p_index] = tuples_recv_right[i];
    t2[p_index] = tuples_recv_right[i+1];
    t3[p_index] = tuples_recv_right[i+2];
  }
}
