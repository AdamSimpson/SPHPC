#ifndef SPH_SRC_GEOMETRY_H_
#define SPH_SRC_GEOMETRY_H_

// Forward Declaration
struct Particles;
struct Params;
struct AABB;
struct Communication;

//Axis aligned bounding box
struct AABB {
  double min_x;
  double max_x;
  double min_y;
  double max_y;
  double min_z;
  double max_z;
};

// Places fluid particles equally spaced within the fluid AABB
void ConstructFluidVolume(struct Particles *const particles,
                          struct Params *const params,
                          const struct AABB *const fluid);

// Sets fluid counts
void SetParticleNumbers(const struct AABB *const fluid_global,
                        struct Particles *particles,
                        struct Params *const params,
                        struct Communication *const communication);

// Sets node x max and x min based upon particle count
void CheckPartition(const struct Particles *particles,
                    struct Params *const params);

#endif
