#include "communication.h"

#include <stdio.h>
#include <stddef.h>
#include <string.h>

void allocate_communication(communication_t *const communication)
{
    // Allocate index arrays
    size_t index_bytes = communication->max_comm_particles*sizeof(int);

    communication->edges.edge_indices_left  = malloc(index_bytes);
    communication->edges.edge_indices_right = malloc(index_bytes);

    communication->out_of_bounds.oob_indices_left = malloc(index_bytes);
    communication->out_of_bounds.oob_indices_right = malloc(index_bytes);

    // Allocate send and receive buffers
    int num_components=3;
    size_t particle_bytes = communication->max_comm_particles*sizeof(fluid_particle_t);
    size_t component_bytes = num_components*communication->max_comm_particles*sizeof(double);

    communication->particle_send_buffer = malloc(particle_bytes);

    communication->halo_components_send_buffer = malloc(component_bytes);
    communication->halo_components_recv_buffer = malloc(component_bytes);
}

void free_communication(communication_t *const communication)
{
    free(communication->edges.edge_indices_left);
    free(communication->edges.edge_indices_right);
    free(communication->out_of_bounds.oob_indices_left);
    free(communication->out_of_bounds.oob_indices_right);
    free(communication->particle_send_buffer);
    free(communication->halo_components_send_buffer);
    free(communication->halo_components_recv_buffer);
}

void init_communication(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    create_MPI_types();
}

void finalize_communication()
{
  free_MPI_types();
  MPI_Finalize();
}

int get_rank()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}

int get_num_procs()
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  return nprocs;
}

void create_MPI_types()
{
    //Create fluid particle type
    MPI_Datatype types[18];
    int i;
    for (i=0; i<17; i++) types[i] = MPI_DOUBLE;
    types[17] = MPI_INT;
    int blocklens[18];
    for (i=0; i<18; i++) blocklens[i] = 1;
    MPI_Aint disps[18];
    // Get displacement of each struct member
    disps[0] = offsetof( fluid_particle_t,  x_star);
    disps[1] = offsetof( fluid_particle_t,  y_star);
    disps[2] = offsetof( fluid_particle_t,  z_star);
    disps[3] = offsetof( fluid_particle_t,  x);
    disps[4] = offsetof( fluid_particle_t,  y);
    disps[5] = offsetof( fluid_particle_t,  z);
    disps[6] = offsetof( fluid_particle_t,  v_x);
    disps[7] = offsetof( fluid_particle_t,  v_y);
    disps[8] = offsetof( fluid_particle_t,  v_z);
    disps[9] = offsetof( fluid_particle_t,  dp_x);
    disps[10] = offsetof( fluid_particle_t, dp_y);
    disps[11] = offsetof( fluid_particle_t, dp_z);
    disps[12] = offsetof( fluid_particle_t, w_x);
    disps[13] = offsetof( fluid_particle_t, w_y);
    disps[14] = offsetof( fluid_particle_t, w_z);
    disps[15] = offsetof( fluid_particle_t, density);
    disps[16] = offsetof( fluid_particle_t, lambda);
    disps[17] = offsetof( fluid_particle_t, id);
    // Commit type
    MPI_Type_create_struct( 18, blocklens, disps, types, &Particletype );
    MPI_Type_commit( &Particletype );
}

void free_MPI_types()
{
    MPI_Type_free(&Particletype);
}

void start_halo_exchange(communication_t *const communication,
                         const param_t *const params,
                         fluid_particle_t *const fluid_particles)
{
    int i;
    int rank = params->rank;
    int nprocs = params->nprocs;
    edge_t *const edges = &communication->edges;

    const fluid_particle_t *p;
    double h = params->smoothing_radius;

    // Set edge particle indices and update number
    edges->number_edge_particles_left = 0;
    edges->number_edge_particles_right = 0;
    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p = &fluid_particles[i];
        if (p->x_star - params->node_start_x <= h)
            edges->edge_indices_left[edges->number_edge_particles_left++] = i;
        else if (params->node_end_x - p->x_star <= h)
            edges->edge_indices_right[edges->number_edge_particles_right++] = i;
    }

    int num_moving_left = edges->number_edge_particles_left;
    int num_moving_right = edges->number_edge_particles_right;

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    // Get number of halo particles from right and left
    int num_from_left = 0;
    int num_from_right = 0;
    int tag = 3217;
    // Send number to right and receive from left
    MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag, &num_from_left,1,MPI_INT,proc_to_left,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    // Send number to left and receive from right
    tag = 8425;
    MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left, tag, &num_from_right,1,MPI_INT,proc_to_right,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    printf("total moving: %d\n", num_moving_right + num_moving_left);

    // Set send buffer points
    fluid_particle_t *sendl_buffer = communication->particle_send_buffer;
    fluid_particle_t *sendr_buffer = sendl_buffer + num_moving_left;

    // Pack send buffers
    for(i=0; i<edges->number_edge_particles_left; i++)
        sendl_buffer[i] = fluid_particles[edges->edge_indices_left[i]];
    for(i=0; i<edges->number_edge_particles_right; i++)
        sendr_buffer[i] = fluid_particles[edges->edge_indices_right[i]];

    //Index to start receiving halo particles
    int indexToReceiveLeft = params->number_fluid_particles_local;
    int indexToReceiveRight = indexToReceiveLeft + num_from_left;

    int tagl = 4312;
    int tagr = 5177;

    // Receive halo from left rank
    MPI_Irecv(&fluid_particles[indexToReceiveLeft], num_from_left, Particletype, proc_to_left,tagl, MPI_COMM_WORLD, &edges->reqs[0]);
    // Receive halo from right rank
    MPI_Irecv(&fluid_particles[indexToReceiveRight], num_from_right, Particletype, proc_to_right,tagr, MPI_COMM_WORLD, &edges->reqs[1]);
    // Send halo to right rank
    MPI_Isend(sendr_buffer, num_moving_right, Particletype, proc_to_right, tagl, MPI_COMM_WORLD, &edges->reqs[2]);
    // Send halo to left rank
    MPI_Isend(sendl_buffer, num_moving_left, Particletype, proc_to_left, tagr, MPI_COMM_WORLD, &edges->reqs[3]);
}

void finish_halo_exchange(communication_t *const communication,
                          const fluid_particle_t *const fluid_particles,
                          param_t *const params)
{
    // Wait for transfer to complete
    MPI_Status statuses[4];
    MPI_Waitall(4, communication->edges.reqs, statuses);

    int num_received_right = 0;
    int num_received_left = 0;
    MPI_Get_count(&statuses[0], Particletype, &num_received_left);
    MPI_Get_count(&statuses[1], Particletype, &num_received_right);

    params->number_halo_particles_left  = num_received_left;
    params->number_halo_particles_right = num_received_right;

    printf("rank %d, halo: recv %d from left, %d from right\n", params->rank,num_received_left,num_received_right);
}

// Transfer particles that are out of node bounds
void transfer_OOB_particles(const communication_t *const communication,
                            fluid_particle_t *const fluid_particles,
                            param_t *const params)
{
    int i;
    const fluid_particle_t *p;
    int rank = params->rank;
    int nprocs = params->nprocs;

    int i_left = 0;
    int i_right = 0;

    // Set send buffers
    fluid_particle_t *const sendl_buffer = communication->particle_send_buffer;
    fluid_particle_t *const sendr_buffer = sendl_buffer + communication->max_comm_particles;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = &fluid_particles[i];

        // Set OOB particle indices and update number
        // remove particle from end of array to fill vacancy
        if (p->x_star < params->node_start_x && params->rank != 0) {
            sendl_buffer[i_left++] = *p;
            fluid_particles[i] = fluid_particles[params->number_fluid_particles_local-1];
            fluid_particles[i].id = i;
            params->number_fluid_particles_local--;
        }
        else if (p->x_star > params->node_end_x && params->rank != params->nprocs-1) {
            sendr_buffer[i_right++] = *p;
            fluid_particles[i] = fluid_particles[params->number_fluid_particles_local-1];
            fluid_particles[i].id = i;
            params->number_fluid_particles_local--;
        }
    }

    int num_moving_left = i_left;
    int num_moving_right = i_right;

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    // Get number of particles from right and left
    int num_from_left = 0;
    int num_from_right = 0;
    int tag = 7006;
    // Send number to right and receive from left
    MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag, &num_from_left,1,MPI_INT,proc_to_left,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    // Send number to left and receive from right
    tag = 8278;
    MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left, tag, &num_from_right,1,MPI_INT,proc_to_right,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Allocate recieve buffers
    fluid_particle_t *recvl_buffer = &fluid_particles[params->number_fluid_particles_local];
    fluid_particle_t *recvr_buffer = &fluid_particles[params->number_fluid_particles_local + num_from_left];

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

    int total_received = num_received_right + num_received_left;

    // Update received particle ID
    for(i=params->number_fluid_particles_local; i<params->number_fluid_particles_local+total_received; i++)
    {
        fluid_particles[i].id = i;
    }

    // Update number of particles
    params->number_fluid_particles_local += total_received;

    printf("rank %d OOB: sent left %d, right: %d recv left:%d, right: %d\n", rank, num_moving_left, num_moving_right, num_received_left, num_received_right);
}

void update_halo_lambdas(const communication_t *const communication,
                         const param_t *const params,
                         fluid_particle_t *const fluid_particles)
{
    int i;
    fluid_particle_t *p;
    const edge_t *const edges = &communication->edges;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int num_moving_left  = edges->number_edge_particles_left;
    int num_moving_right = edges->number_edge_particles_right;

    int num_from_left  = params->number_halo_particles_left;
    int num_from_right = params->number_halo_particles_right;

    // Allocate send/recv buffers
    double *send_lambdas_left  = communication->halo_components_send_buffer;
    double *send_lambdas_right = send_lambdas_left + num_moving_left;

    double *recv_lambdas_left  = communication->halo_components_recv_buffer;
    double *recv_lambdas_right = recv_lambdas_left + num_from_left;

    // Pack local halo lambdas
    for(i=0; i<num_moving_left; i++) {
        p = &fluid_particles[edges->edge_indices_left[i]];
        send_lambdas_left[i] = p->lambda;
    }
    for(i=0; i<num_moving_right; i++) {
        p = &fluid_particles[edges->edge_indices_right[i]];
        send_lambdas_right[i] = p->lambda;
    }

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

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
    for(i=0; i<num_from_left; i++) {
        p = &fluid_particles[params->number_fluid_particles_local + i];;
        p->lambda = recv_lambdas_left[i];
    }
    for(i=0; i<num_from_right; i++) {
        p = &fluid_particles[ params->number_fluid_particles_local + num_from_left + i];
        p->lambda = recv_lambdas_right[i];
    }
}

void update_halo_positions(const communication_t *const communication,
                           const param_t *const params,
                           fluid_particle_t *const fluid_particles)
{
    int i;
    fluid_particle_t *p;
    const edge_t *const edges = &communication->edges;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // x,y,z components required
    int num_components = 3;
    int num_moving_left = num_components*edges->number_edge_particles_left;
    int num_moving_right = num_components*edges->number_edge_particles_right;

    int num_from_left = num_components*params->number_halo_particles_left;
    int num_from_right = num_components*params->number_halo_particles_right;

    // Allocate send/recv buffers
    // Could combine left/right into single malloc...
    double *send_positions_left = communication->halo_components_send_buffer;
    double *send_positions_right = send_positions_left + num_moving_left;

    double *recv_positions_left = communication->halo_components_recv_buffer;
    double *recv_positions_right = recv_positions_left + num_from_left;

    // Pack local edge positions
    for(i=0; i<num_moving_left; i+=3) {
        p = &fluid_particles[edges->edge_indices_left[i/3]];
        send_positions_left[i]   = p->x_star;
        send_positions_left[i+1] = p->y_star;
        send_positions_left[i+2] = p->z_star;
    }
    for(i=0; i<num_moving_right; i+=3) {
        p = &fluid_particles[edges->edge_indices_right[i/3]];
        send_positions_right[i]   = p->x_star;
        send_positions_right[i+1] = p->y_star;
        send_positions_right[i+2] = p->z_star;
    }

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    // Could do async to perhaps increase performance
    // Send positions to right and receive from left
    int tag = 874;
    MPI_Sendrecv(send_positions_right, num_moving_right, MPI_DOUBLE, proc_to_right, tag,
                 recv_positions_left, num_from_left, MPI_DOUBLE, proc_to_left, tag,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Send positions to left and receive from right
    tag = 546;
    MPI_Sendrecv(send_positions_left, num_moving_left, MPI_DOUBLE, proc_to_left, tag,
                 recv_positions_right, num_from_right, MPI_DOUBLE, proc_to_right,tag,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Unpack halo particle positions
    for(i=0; i<num_from_left; i+=3) {
        p = &fluid_particles[params->number_fluid_particles_local + i/3];;
        p->x_star = recv_positions_left[i];
        p->y_star = recv_positions_left[i+1];
        p->z_star = recv_positions_left[i+2];
    }
    for(i=0; i<num_from_right; i+=3) {
        p = &fluid_particles[ params->number_fluid_particles_local + num_from_left/3 + i/3];
        p->x_star = recv_positions_right[i];
        p->y_star = recv_positions_right[i+1];
        p->z_star = recv_positions_right[i+2];
    }
}
