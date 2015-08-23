#define restrict __restrict

#include "obstacle.h"
#include "print_macros.h"
#include "cuda_runtime.h"
#include "particles.h"
#include "obstacle.h"
#include "obstacle_cuda.h"

inline __device__ float dot(float3 a, float3 b) {
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline __device__ float3 normalize(float3 v) {
    float inv_length = rsqrtf(dot(v, v));
    return make_float3(v.x*inv_length, v.y*inv_length, v.z*inv_length);
}


// Test if inside obstacle bounding box
__device__ bool InsideObstacleBounds(const struct Obstacle obstacle,
                               const double x, const double y, const double z) {

    struct AABB world_bounds = obstacle.world_bounds;

    if(x > world_bounds.min_x && x < world_bounds.max_x &&
       y > world_bounds.min_y && y < world_bounds.max_y &&
       z > world_bounds.min_z && z < world_bounds.max_z)
        return true;
    else
      return false;
}

inline __device__ float3 CalculateNormal(cudaTextureObject_t SDF_tex,
                                  const float x_tex,
                                  const float y_tex,
                                  const float z_tex,
                                  const float x_spacing,
                                  const float y_spacing,
                                  const float z_spacing) {
  // Calculate gradient of SDF(surface normal)
  const float x_right = tex3D<float>(SDF_tex, x_tex + x_spacing, y_tex, z_tex);
  const float x_left  = tex3D<float>(SDF_tex, x_tex - x_spacing, y_tex, z_tex);
  const float y_up    = tex3D<float>(SDF_tex, x_tex, y_tex + y_spacing, z_tex);
  const float y_down  = tex3D<float>(SDF_tex, x_tex, y_tex - y_spacing, z_tex);
  const float z_front = tex3D<float>(SDF_tex, x_tex, y_tex, z_tex + z_spacing);
  const float z_back  = tex3D<float>(SDF_tex, x_tex, y_tex, z_tex - z_spacing);

  float3 gradient = make_float3(x_right - x_left,
                                y_up - y_down,
                                z_front - z_back);

  return normalize(gradient);
}

// Test if particle is inside obstacle and if so remove it
__global__ void CalculateParticleCollisions_kernel(double *restrict x,
                                                   double *restrict y,
                                                   double *restrict z,
                                                   const int particle_count,
                                                   const struct Obstacle obstacle) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if(i < particle_count) {
    // World coordinates
    const double x_world = x[i];
    const double y_world = y[i];
    const double z_world = z[i];

    if(InsideObstacleBounds(obstacle, x_world, y_world, z_world)) {

      const struct AABB obstacle_bounds_world = obstacle.world_bounds;

      // Compute normalized [0,1] tex coordinates
      const float obstacle_length_world = (obstacle_bounds_world.max_x
                                          - obstacle_bounds_world.min_x);
      const float x_tex = (x_world - obstacle_bounds_world.min_x)
                          /obstacle_length_world;
      const float obstacle_height_world = (obstacle_bounds_world.max_y
                                          - obstacle_bounds_world.min_y);
      const float y_tex = (y_world - obstacle_bounds_world.min_y)
                          /obstacle_height_world;
      const float obstacle_depth_world = (obstacle_bounds_world.max_z
                                         - obstacle_bounds_world.min_z);
      const float z_tex = (z_world - obstacle_bounds_world.min_z)
                          /obstacle_depth_world;

      // Fetch model space distance to surface from texture
      const float surface_distance_model = tex3D<float>(obstacle.obstacle_cuda.SDF_tex, x_tex, y_tex, z_tex);

      if(surface_distance_model < 0) {  // If inside object
        // SDF grid spacing in normalized texture coordinates [0,1]
        const float x_spacing = 1.0 / (float)obstacle.SDF_dim_x;
        const float y_spacing = 1.0 / (float)obstacle.SDF_dim_y;
        const float z_spacing = 1.0 / (float)obstacle.SDF_dim_z;

        const float3 surface_normal_model = CalculateNormal(obstacle.obstacle_cuda.SDF_tex,
                                                            x_tex, y_tex, z_tex,
                                                            x_spacing, y_spacing, z_spacing);

        // This may not be the exact model length as it is an integer multiple of SDF_spacing
        const float obstacle_length_model = obstacle.SDF_dim_x * obstacle.SDF_spacing;

        // Convert distance and normals in model space to world space
        // This assumes aspect ratio is the same between model and world
        const float model_to_world = obstacle_length_world / obstacle_length_model;
        const float3 surface_normal_world = surface_normal_model;
        const float surface_distance_world = fabsf(surface_distance_model)
                                           * model_to_world;

        y[i] += surface_distance_world * surface_normal_world.y;
        x[i] += surface_distance_world * surface_normal_world.x;
        z[i] += surface_distance_world * surface_normal_world.z;
      } // Inside obstacle
    } // Inside obstacle bounding box
  }
}

extern "C" void CalculateParticleCollisions(double *restrict x,
                                            double *restrict y,
                                            double *restrict z,
                                            const int particle_count,
                                            const struct Obstacle obstacle) {
  int blockSize = 1024;
  int gridSize = (int)ceil((float)particle_count/blockSize);

  CalculateParticleCollisions_kernel<<<gridSize, blockSize>>>(x, y, z,
                                                              particle_count,
                                                              obstacle);

}

void AllocInitObstacle_CUDA(struct Obstacle *obstacle) {

  struct Obstacle_CUDA* obstacle_cuda = &(obstacle->obstacle_cuda);

  // Allocate 3D array for SDF
  const cudaExtent SDF_size = make_cudaExtent(obstacle->SDF_dim_x,
                                              obstacle->SDF_dim_y,
                                              obstacle->SDF_dim_z);

  // Single channel of single precision float format
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
  tex_descriptor.normalizedCoords = 1;
  cudaCreateTextureObject(&(obstacle_cuda->SDF_tex), &resource_descriptor,
                          &tex_descriptor, NULL);

}

void FinalizeObstacle_CUDA(struct Obstacle_CUDA *obstacle_cuda) {
  cudaFreeArray(obstacle_cuda->SDF_buffer);
  cudaDestroyTextureObject(obstacle_cuda->SDF_tex);
}
