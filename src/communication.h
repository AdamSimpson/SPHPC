#ifndef SPH_SRC_COMMUNICATION_H_
#define SPH_SRC_COMMUNICATION_H_

typedef struct EDGES edge_t;
typedef struct OOB oob_t;
typedef struct COMMUNICATION communication_t;

#include "fluid.h"
#include "mpi.h"
#include "simulation.h"

MPI_Datatype Particletype;

// Particles that are within halo width of node edge
struct EDGES {
    int *edge_indices_left;
    int *edge_indices_right;
    int number_edge_particles_left;
    int number_edge_particles_right;
    MPI_Request reqs[4];
};

// Particles that have left the node
struct OOB {
    int *oob_indices_left; // Indicies in particle pointer array for particles traveling left
    int *oob_indices_right;
    int number_oob_particles_left;
    int number_oob_particles_right;
};

struct COMMUNICATION {
    edge_t edges;
    oob_t out_of_bounds;
    int max_comm_particles;
    fluid_particle_t *particle_send_buffer;
    double *halo_components_send_buffer;
    double *halo_components_recv_buffer;
};

void init_communication(int argc, char *argv[]);
void finalize_communication();
void allocate_communication(communication_t *const communication);
void free_communication(communication_t *const communication);
int get_rank();
int get_num_procs();
void create_MPI_types();
void free_MPI_types();
void transfer_OOB_particles(const communication_t *const communication,
                            fluid_particle_t *const fluid_particles,
                            param_t *const params);
void start_halo_exchange(communication_t *const communication,
                         const param_t *const params,
                         fluid_particle_t *const fluid_particles);
void finish_halo_exchange(communication_t *const communication,
                          const fluid_particle_t *const fluid_particles,
                          param_t *const params);
void update_halo_lambdas(const communication_t *const communication,
                         const param_t *const params,
                        fluid_particle_t *const fluid_particles);
void update_halo_positions(const communication_t *const communication,
                          const param_t *const params,
                          fluid_particle_t *const fluid_particles);
#endif
