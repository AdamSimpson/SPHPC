#ifndef SPH_SRC_INPUT_PARSER_H_
#define SPH_SRC_INPUT_PARSER_H_

// Forward Declaration
struct Params;
struct Particles;
struct AABB;
struct Obstacle;

// These functions are to be visible to both C and C++
#ifdef __cplusplus
extern "C" {
#endif

// Read an external parameters file and set simulation values
void ReadParameters(struct Params *parameters,
                    struct Particles *particles,
                    struct AABB *boundary_volume,
                    struct AABB *initial_fluid_volume,
                    struct Obstacle *obstacle);

#ifdef __cplusplus
}
#endif

#endif
