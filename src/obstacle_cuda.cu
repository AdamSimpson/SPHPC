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
  int i = blockIdx.x *blockDim.x + threadIdx.x;
  if(i < particle_count) {
    // World coordinates
    const double world_x = x[i];
    const double world_y = y[i];
    const double world_z = z[i];

    if(InsideObstacleBounds(obstacle, world_x, world_y, world_z)) {

      const struct AABB obstacle_world_bounds = obstacle.world_bounds;

      // Compute normalized tex coordinates
      const float obstacle_world_length = (obstacle_world_bounds.max_x
                                          - obstacle_world_bounds.min_x);
      const float x_tex = (world_x - obstacle_world_bounds.min_x)/obstacle_world_length;

      const float obstacle_world_height = (obstacle_world_bounds.max_y
                                          - obstacle_world_bounds.min_y);
      const float y_tex = (world_y - obstacle_world_bounds.min_y)/obstacle_world_height;

      const float obstacle_world_depth = (obstacle_world_bounds.max_z
                                         - obstacle_world_bounds.min_z);
      const float z_tex = (world_z - obstacle_world_bounds.min_z)/obstacle_world_depth;

      // Fetch model space distance to surface from texture
      const float surface_model_distance = tex3D<float>(obstacle.obstacle_cuda.SDF_tex, x_tex, y_tex, z_tex);

      if(surface_model_distance < 0) {  // If inside object
        // SDF grid spacing in normalized texture coordinates [0,1]
        const float x_spacing = 1.0 / (float)obstacle.SDF_dim_x;
        const float y_spacing = 1.0 / (float)obstacle.SDF_dim_y;
        const float z_spacing = 1.0 / (float)obstacle.SDF_dim_z;

        const float3 model_surface_normal = CalculateNormal(obstacle.obstacle_cuda.SDF_tex,
                                                            x_tex, y_tex, z_tex,
                                                            x_spacing, y_spacing, z_spacing);

        // This may not be the exact model length as it is an integer multiple of SDF_spacing
        const float obstacle_model_length = obstacle.SDF_dim_x * obstacle.SDF_spacing;

        // Convert distance and normals in model space to world space
        // This assumes aspect ratio is the same between model and world
        const float model_to_world = obstacle_world_length / obstacle_model_length;
        const float3 surface_world_normal = model_surface_normal;
        const float surface_world_distance = fabsf(surface_model_distance) * model_to_world;

        x[i] += surface_world_distance * surface_world_normal.x;
        y[i] += surface_world_distance * surface_world_normal.y;
        z[i] += surface_world_distance * surface_world_normal.z;
      }

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

  CalculateParticleCollisions_kernel<<<gridSize, blockSize>>>(x,y,z,
                                                              particle_count,
                                                              obstacle);

}

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
  tex_descriptor.normalizedCoords = 1;
  cudaCreateTextureObject(&(obstacle_cuda->SDF_tex), &resource_descriptor,
                          &tex_descriptor, NULL);

}

void FinalizeObstacle_CUDA(struct Obstacle_CUDA *obstacle_cuda) {
  cudaFreeArray(obstacle_cuda->SDF_buffer);
  cudaDestroyTextureObject(obstacle_cuda->SDF_tex);
}
