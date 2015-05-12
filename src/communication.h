#ifndef SPH_SRC_COMMUNICATION_H_
#define SPH_SRC_COMMUNICATION_H_

typedef struct EDGES Edges;
typedef struct OOB OOB;
typedef struct COMMUNICATION Communication;

#include "mpi.h"
#include "fluid.h"
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
  Edges edges;
  OOB out_of_bounds;
  int max_comm_particles;
  FluidParticle *particle_send_buffer;
  double *halo_components_send_buffer;
  double *halo_components_recv_buffer;
};

void InitCommunication(int argc, char *argv[]);
void FinalizeCommunication();
void AllocateCommunication(Communication *const communication);
void FreeCommunication(Communication *const communication);
int get_rank();
int get_num_procs();
void CreateMPITypes();
void FreeMPITypes();

void TransferOOBParticles(const Communication *const communication,
                          FluidParticle *const fluid_particles,
                          Params *const params);
void StartHaloExchange(Communication *const communication,
                       const Params *const params,
                       FluidParticle *const fluid_particles);
void FinishHaloExchange(Communication *const communication,
                        const FluidParticle *const fluid_particles,
                        Params *const params);
void UpdateHaloLambdas(const Communication *const communication,
                       const Params *const params,
                       FluidParticle *const fluid_particles);
void UpdateHaloPositions(const Communication *const communication,
                         const Params *const params,
                         FluidParticle *const fluid_particles);
#endif
