#ifndef SPH_SRC_GEOMETRY_H_
#define SPH_SRC_GEOMETRY_H_

// Forward Declaration
struct FluidParticles;
struct Params;
struct AABB;
struct Communication;

struct AABB {
  double min_x;
  double max_x;
  double min_y;
  double max_y;
  double min_z;
  double max_z;
};//Axis aligned bounding box

void ConstructFluidVolume(struct FluidParticles *const particles,
                          struct Params *const params,
                          const struct AABB *const fluid);

void SetParticleNumbers(const struct AABB *const fluid_global,
                        struct Params *const params,
                        struct Communication *const communication);

void CheckPartition(struct Params *const params);

#endif
