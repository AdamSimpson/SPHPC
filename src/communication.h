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
  int *indices_left;
  int *indices_right;
  int number_particles_left;
  int number_particles_right;
};

// Particles that have left the node
struct OOB {
  int *indices_left; // Indicies in particle pointer array for particles traveling left
  int *indices_right;
  int number_particles_left;
  int number_particles_right;
};

struct COMMUNICATION {
  Edges edges;
  OOB out_of_bounds;
  int max_comm_particles;
  double *particle_send_buffer;
  double *particle_recv_buffer;
};

void InitCommunication(int argc, char *argv[]);
void FinalizeCommunication();
void AllocateCommunication(Communication *const communication);
void FreeCommunication(Communication *const communication);
int get_rank();
int get_num_procs();
void PackHaloComponents(const Communication *const communication,
                        const FluidParticles *const fluid_particles,
                        double *const packed_send_left,
                        double *const packed_send_right);

void UnpackHaloComponents(const Params *const params,
                          const double *const packed_recv_left,
                          const double *const packed_recv_right,
                          FluidParticles *const fluid_particles);

void PackOOBComponents(const Communication *const communication,
                       const FluidParticles *const fluid_particles,
                       double *const packed_send_left,
                       double *const packed_send_right);

void UnpackOOBComponents(const int num_from_left, const int num_from_right,
                         const double *const packed_recv_left,
                         const double *const packed_recv_right,
                         Params *const params,
                         FluidParticles *const fluid_particles);

void TransferOOBParticles(Communication *const communication,
                          FluidParticles *const fluid_particles,
                          Params *const params);

void HaloExchange(Communication *const communication,
                  Params *const params,
                  FluidParticles *const fluid_particles);

void UpdateHaloLambdas(const Communication *const communication,
                       const Params *const params,
                       FluidParticles *const fluid_particles);
void UpdateHaloPositions(const Communication *const communication,
                         const Params *const params,
                         FluidParticles *const fluid_particles);
#endif
