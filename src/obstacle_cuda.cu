#include "obstacle.h"
#include "print_macros.h"
#include "cuda_runtime.h"
#include "obstacle_cuda.h"

void AllocInitObstacle_CUDA(struct Obstacle *obstacle) {

  struct Obstacle_CUDA* obstacle_cuda = &(obstacle->obstacle_cuda);

  // Allocate 3D array for SDF
  const cudaExtent SDF_size = make_cudaExtent(obstacle->SDF_dim_x,
                                              obstacle->SDF_dim_y,
                                              obstacle->SDF_dim_z);

  const cudaChannelFormatDesc channel_desc = cudaCreateChannelDesc(
                                               32,
                                               0, 0, 0,
                                               cudaChannelFormatKindFloat);
  cudaError_t err;
  err = cudaMalloc3DArray(&(obstacle_cuda->SDF_buffer), &channel_desc, SDF_size);
  if(err != cudaSuccess)
    EXIT_PRINT("cudaMalloc3DArray error: %d!\n", err);

  // Copy 3D array to GPU
  cudaMemcpy3DParms copy_params = {0};
  copy_params.srcPtr   = make_cudaPitchedPtr((void*)obstacle->SDF_buffer,
                                             SDF_size.width*sizeof(float),
                                             SDF_size.width,
                                             SDF_size.height);
  copy_params.dstArray = obstacle_cuda->SDF_buffer;
  copy_params.extent   = SDF_size;
  copy_params.kind     = cudaMemcpyHostToDevice;
  err = cudaMemcpy3D(&copy_params);
  if(err != cudaSuccess)
    EXIT_PRINT("cudaMalloc3DArray error! %d\n", err);

  // Create texture object
  cudaResourceDesc resource_descriptor;
  memset(&resource_descriptor, 0, sizeof(resource_descriptor));
  resource_descriptor.resType = cudaResourceTypeArray;
  resource_descriptor.res.array.array = obstacle_cuda->SDF_buffer;

  // Create texture object
  cudaTextureDesc tex_descriptor;
  memset(&tex_descriptor, 0, sizeof(tex_descriptor));
  tex_descriptor.filterMode = cudaFilterModeLinear;
  tex_descriptor.addressMode[0] = cudaAddressModeClamp;
  tex_descriptor.addressMode[1] = cudaAddressModeClamp;
  tex_descriptor.addressMode[2] = cudaAddressModeClamp;
  cudaCreateTextureObject(&(obstacle_cuda->SDF_tex), &resource_descriptor,
                          &tex_descriptor, NULL);

}

void FinalizeObstacle_CUDA(struct Obstacle_CUDA *obstacle_cuda) {
  cudaFreeArray(obstacle_cuda->SDF_buffer);
  cudaDestroyTextureObject(obstacle_cuda->SDF_tex);
}
/*
// Test if inside obstacle bounding box
bool InsideObstacle(const struct Obstacle *obstacle,
                    const double x, const double y, const double z) {

    struct AABB world_bounds = obstacle->world_bounds;
    if(x < world_bounds.min_x || x > world_bounds.max_x)
        return true;
    if(y < world_bounds.min_y || y > world_bounds.max_y)
        return true;
    if(z < world_bounds.min_z || z > world_bounds.max_z)
        return true;
}
*/
