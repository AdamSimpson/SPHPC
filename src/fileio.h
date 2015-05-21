#ifndef SPH_SRC_FILEIO_H_
#define SPH_SRC_FILEIO_H_

// Forward Declaration
struct Particles;
struct Params;

struct FileIO {
    double *write_buffer;
    char *output_path;
    int file_num;
};

// Allocate output file buffers
void FileIOInit(struct FileIO *const file_io,
                const struct Particles *const particles,
                const struct Params *params);

// Free file buffers
void FileIOFinalize(struct FileIO *const file_io);

// Writes particle x,y,z positions to file using MPI-IO
void WriteMPI(const struct Particles *const particles,
              const struct Params *const params,
              struct FileIO *const file_io);

#endif
