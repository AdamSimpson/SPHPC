#ifndef SPH_SRC_FILEIO_H_
#define SPH_SRC_FILEIO_H_

#include "fluid.h"
#include "simulation.h"

void WriteMPI(const FluidParticles *const particles,
              const Params *const params,
              const int fileNum);

#endif
