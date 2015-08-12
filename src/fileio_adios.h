#ifndef SPH_SRC_FILEIO_ADIOS_H_
#define SPH_SRC_FILEIO_ADIOS_H_

#include <inttypes.h>

// Forward Declaration
struct Particles;
struct Params;

struct FileIO {
    char *output_path;
    int file_num;
    int rank;
    int64_t adios_handle;
};

// Allocate output file buffers
void AllocInitFileIO(struct FileIO *const file_io,
                const struct Particles *const particles,
                const struct Params *params);

// Free file buffers
void FinalizeFileIO(struct FileIO *const file_io);

// Writes particle x,y,z positions to file using MPI-IO
void WriteParticles(const struct Particles *const particles,
              const struct Params *const params,
              struct FileIO *const file_io);

#endif
