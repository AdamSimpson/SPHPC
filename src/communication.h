#ifndef fluid_communication_h
#define fluid_communication_h

typedef struct EDGES edge_t;
typedef struct OOB oob_t;
typedef struct COMMUNICATION communication_t;

#include "fluid.h"
#include "mpi.h"
#include "simulation.h"

MPI_Datatype Particletype;

// Particles that are within halo width of node edge
struct EDGES {
    int max_edge_particles;
    int *edge_indices_left;
    int *edge_indices_right;
    int number_edge_particles_left;
    int number_edge_particles_right;

    MPI_Request reqs[4];
};

// Particles that have left the node
struct OOB {
    int max_oob_particles;
    int *oob_indices_left; // Indicies in particle pointer array for particles traveling left
    int *oob_indices_right;
    int number_oob_particles_left;
    int number_oob_particles_right;
};

struct COMMUNICATION {
    edge_t edges;
    oob_t out_of_bounds;
    fluid_particle_t *particle_send_buffer;
    double *halo_components_send_buffer;
    double *halo_components_recv_buffer;
};

void init_communication(int argc, char *argv[]);
void allocate_communication(communication_t *communication);
int get_rank();
int get_num_procs();
void createMpiTypes();
void freeMpiTypes();
void transferHalos(communication_t *communication, fluid_particle_t *fluid_particles, param_t *params);
void transferOOBParticles(communication_t *communication, fluid_particle_t *fluid_particles, param_t *params);
void startHaloExchange(communication_t *communication, fluid_particle_t *fluid_particles, param_t *params);
void finishHaloExchange(communication_t *communication, fluid_particle_t *fluid_particles, param_t *params);
void update_halo_lambdas(communication_t *communication, fluid_particle_t *fluid_particles, param_t *params);
void update_halo_positions(communication_t *communication, fluid_particle_t *fluid_particles, param_t *params);

#endif
