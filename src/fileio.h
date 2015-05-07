#ifndef SPH_SRC_FILEIO_H_
#define SPH_SRC_FILEIO_H_

#include "fluid.h"
#include "simulation.h"

void write_MPI(fluid_particle_t *particles, int fileNum, param_t *params);

#endif
