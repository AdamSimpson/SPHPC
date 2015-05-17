#ifndef SPH_SRC_FILEIO_H_
#define SPH_SRC_FILEIO_H_

// Forward Declaration
struct FluidParticles;
struct Params;

void WriteMPI(const struct FluidParticles *const particles,
              const struct Params *const params,
              const int fileNum);

#endif
