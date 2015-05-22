#ifndef SPH_SRC_COMMUNICATION_H_
#define SPH_SRC_COMMUNICATION_H_

// Forward Declaration
struct Particles;
struct Params;

// Contains indices to particles that are within smooth radius of node edge
struct Edges {
  int *restrict indices_left;
  int *restrict indices_right;
  int particle_count_left;
  int particle_count_right;
};

// Contains indicies to particles that have gone Out Of Bounds
struct OOB {
  int *restrict indices_left;
  int *restrict indices_right;
  int particle_count_left;
  int particle_count_right;
};

// Contains MPI send/recv buffers and above structs
struct Communication {
  struct Edges edges;
  struct OOB out_of_bounds;
  int max_particles;
  double *restrict particle_send_buffer;
  double *restrict particle_recv_buffer;
};

// Wrapper for MPI_Init
void InitCommunication(int argc, char *argv[]);

// Wrapper for MPI_Finalize
void FinalizeCommunication();

// Allocates edge/oob index arrays and MPI send/recv buffers
void AllocateCommunication(struct Communication *const communication);

// Free memory allocated in AllocateCommunications
void FreeCommunication(struct Communication *const communication);

// Return MPI Rank
int get_rank();

// Return MPI processor count
int get_proc_count();

// Packs particles  arraydouble components into packed_send_left
// and packed_send_right based upon indices held in edges
// struct member of communication
void PackHaloComponents(const struct Communication *const communication,
                        const struct Particles *const particles,
                        double *const packed_send_left,
                        double *const packed_send_right);

// appends packed_recv_left followed by packed_recv_right doubles on to
// particles array
void UnpackHaloComponents(const double *const packed_recv_left,
                          const double *const packed_recv_right,
                          struct Particles *const particles);

// Packs particles  arraydouble components into packed_send_left
// and packed_send_right based upon indices held in edges
// struct member of communication
void PackOOBComponents(const struct Communication *const communication,
                       const struct Particles *const particles,
                       double *const packed_send_left,
                       double *const packed_send_right);

// appends packed_recv_left followed by packed_recv_right doubles on to
// particles array
void UnpackOOBComponents(const int num_from_left, const int num_from_right,
                         const double *const packed_recv_left,
                         const double *const packed_recv_right,
                         struct Particles *const particles);

// Identifies particles that have left the node start/end x boundary
// then packs them into flat double arrays and sends them to the
// left or right node
void ExchangeOOB(struct Communication *const communication,
                          struct Particles *const particles,
                          struct Params *const params);

// Identifies particles that have are within one smooth radius of node edge
// then packs them into flat double arrays and sends them to the
// left or right node
void ExchangeHalo(struct Communication *const communication,
                  struct Params *const params,
                  struct Particles *const particles);

// Updates each halo particles scalar lambda value
void UpdateHaloLambdas(const struct Communication *const communication,
                       const struct Params *const params,
                       struct Particles *const particles);

// Updates each halo particles 3 component position value
void UpdateHaloPositions(const struct Communication *const communication,
                         const struct Params *const params,
                         struct Particles *const particles);
#endif
