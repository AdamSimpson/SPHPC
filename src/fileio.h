#ifndef SPH_SRC_FILEIO_H_
#define SPH_SRC_FILEIO_H_

#include "fluid.h"
#include "simulation.h"

void write_MPI(const fluid_particle_t *const particles,
               const param_t *const params,
               const int fileNum);

#endif
