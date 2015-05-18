#ifndef SPH_SRC_COMMUNICATION_H_
#define SPH_SRC_COMMUNICATION_H_

// Forward Declaration
struct FluidParticles;
struct Params;

// Particles that are within halo width of node edge
struct Edges {
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

struct Communication {
  struct Edges edges;
  struct OOB out_of_bounds;
  int max_comm_particles;
  double *particle_send_buffer;
  double *particle_recv_buffer;
};

void InitCommunication(int argc, char *argv[]);
void FinalizeCommunication();
void AllocateCommunication(struct Communication *const communication);
void FreeCommunication(struct Communication *const communication);
int get_rank();
int get_num_procs();
void PackHaloComponents(const struct Communication *const communication,
                        const struct FluidParticles *const particles,
                        double *const packed_send_left,
                        double *const packed_send_right);

void UnpackHaloComponents(const struct Params *const params,
                          const double *const packed_recv_left,
                          const double *const packed_recv_right,
                          struct FluidParticles *const particles);

void PackOOBComponents(const struct Communication *const communication,
                       const struct FluidParticles *const particles,
                       double *const packed_send_left,
                       double *const packed_send_right);

void UnpackOOBComponents(const int num_from_left, const int num_from_right,
                         const double *const packed_recv_left,
                         const double *const packed_recv_right,
                         struct Params *const params,
                         struct FluidParticles *const particles);

void TransferOOBParticles(struct Communication *const communication,
                          struct FluidParticles *const particles,
                          struct Params *const params);

void HaloExchange(struct Communication *const communication,
                  struct Params *const params,
                  struct FluidParticles *const particles);

void UpdateHaloLambdas(const struct Communication *const communication,
                       const struct Params *const params,
                       struct FluidParticles *const particles);
void UpdateHaloPositions(const struct Communication *const communication,
                         const struct Params *const params,
                         struct FluidParticles *const particles);
#endif
