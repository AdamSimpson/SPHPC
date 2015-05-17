#ifndef SPH_SRC_INPUT_PARSER_H_
#define SPH_SRC_INPUT_PARSER_H_

// Forward Declaration
struct Params;
struct AABB;

// These functions are to be visible to both C and C++
#ifdef __cplusplus
extern "C" {
#endif

void ReadParameters(struct Params *const parameters,
                    struct AABB *boundary_volume,
                    struct AABB *initial_fluid_volume);
#ifdef __cplusplus
}
#endif

#endif
