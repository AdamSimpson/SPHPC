#ifndef SPH_SRC_COMMUNICATION_H_
#define SPH_SRC_COMMUNICATION_H_

// Forward Declaration
struct Particles;
struct Params;

// Particles that are within halo width of node edge
struct Edges {
  int *restrict indices_left;
  int *restrict indices_right;
  int particle_count_left;
  int particle_count_right;
};

// Particles that have left the node
struct OOB {
  int *restrict indices_left; // Indicies in particle pointer array for particles traveling left
  int *restrict indices_right;
  int particle_count_left;
  int particle_count_right;
};

struct Communication {
  struct Edges edges;
  struct OOB out_of_bounds;
  int max_particles;
  double *restrict particle_send_buffer;
  double *restrict particle_recv_buffer;
};

void InitCommunication(int argc, char *argv[]);
void FinalizeCommunication();
void AllocateCommunication(struct Communication *const communication);
void FreeCommunication(struct Communication *const communication);
int get_rank();
int get_proc_count();
void PackHaloComponents(const struct Communication *const communication,
                        const struct Particles *const particles,
                        double *const packed_send_left,
                        double *const packed_send_right);

void UnpackHaloComponents(const double *const packed_recv_left,
                          const double *const packed_recv_right,
                          struct Particles *const particles);

void PackOOBComponents(const struct Communication *const communication,
                       const struct Particles *const particles,
                       double *const packed_send_left,
                       double *const packed_send_right);

void UnpackOOBComponents(const int num_from_left, const int num_from_right,
                         const double *const packed_recv_left,
                         const double *const packed_recv_right,
                         struct Particles *const particles);

void TransferOOBParticles(struct Communication *const communication,
                          struct Particles *const particles,
                          struct Params *const params);

void HaloExchange(struct Communication *const communication,
                  struct Params *const params,
                  struct Particles *const particles);

void UpdateHaloLambdas(const struct Communication *const communication,
                       const struct Params *const params,
                       struct Particles *const particles);
void UpdateHaloPositions(const struct Communication *const communication,
                         const struct Params *const params,
                         struct Particles *const particles);
#endif
