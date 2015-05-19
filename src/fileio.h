#ifndef SPH_SRC_FILEIO_H_
#define SPH_SRC_FILEIO_H_

// Forward Declaration
struct FluidParticles;
struct Params;

struct FileIO {
    double *write_buffer;
    char *output_path;
    int file_num;
};

void FileIOInit(struct FileIO *const file_io,
                const struct Params *const params);

void FileIOFinalize(struct FileIO *const file_io);

void WriteMPI(const struct FluidParticles *const particles,
              const struct Params *const params,
              struct FileIO *const file_io);

#endif
