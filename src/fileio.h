#ifndef fluid_fileio_h
#define fluid_fileio_h

#include "fluid.h"
#include "simulation.h"

void write_MPI(fluid_particle_t *particles, int fileNum, param_t *params);

#endif
