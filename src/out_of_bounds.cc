#include "out_of_bounds.h"

#include <thrust/execution_policy.h>
extern "C" {
  #define restrict
  #include "simulation.h"
  #include "communication.h"
  #include "particles.h"
  #undef restrict
}
#include <iostream>

struct LessThan {
  LessThan(double min): x_min(min) {}

  double x_min;

  __host__ __device__
  bool operator()(const double x) {
    return (x < x_min);
  }
};

struct GreaterThan {
  GreaterThan(double max): x_max(max) {}

  double x_max;

  __host__ __device__
  bool operator()(const double x) {
    return (x > x_max);
  }
};

struct OutsideBounds {
  OutsideBounds(double min, double max): x_min(min), x_max(max) {}

  double x_min, x_max;

  __host__ __device__
  bool operator()(const double x) {
    return (x < x_min || x > x_max);
  }
};

// C wrapper function for thrust copy_if with LessThan predicate
extern "C" void CopyIfLessThan(const double min,
                               const int *const input,
                               const int input_count,
                               const double *const stencil,
                               int *const output,
                               int *const output_count) {

  int *end_pointer = thrust::copy_if(thrust::host,
                                     input,
                                     input + input_count,
                                     stencil,
                                     output,
                                     LessThan(min));
  *output_count = end_pointer - output;
}

// C wrapper function for thrust copy_if with GreaterThan predicate
extern "C" void CopyIfGreaterThan(const double max,
                                  const int *const input,
                                  const int input_count,
                                  const double *const stencil,
                                  int *const output,
                                  int *const output_count) {

  int *end_pointer = thrust::copy_if(thrust::host,
                                     input,
                                     input + input_count,
                                     stencil,
                                     output,
                                     GreaterThan(max));
  *output_count = end_pointer - output;
}

// C wrapper function for thrust remove_if with OutsideBounds predicate
extern "C" void RemoveIfOutsideBounds(const double min, const double max,
                                      int *const input,
                                      const int input_count,
                                      const double *const stencil) {

  thrust::remove_if(thrust::host,
                    input,
                    input + input_count,
                    stencil,
                    OutsideBounds(min, max));
}

// Remove particles less than the node start and save to packed array;
void PackOOB(struct Params *const params,
             struct Particles *const particles,
             struct Communication *const communication) {

  // Use x array as stencil, for predicate test, on particle ID array
  // If x is less than node start copy the id to indices_left
  int *end_pointer = thrust::copy_if(thrust::host,
                                     particles->id,
                                     particles->id+particles->local_count,
                                     particles->x,
                                     communication->out_of_bounds.indices_left,
                                     LessThan(params->node_start_x));
  int oob_left_count = end_pointer - communication->out_of_bounds.indices_left;
  communication->out_of_bounds.particle_count_left = oob_left_count;

  // Use x array as stencil, for predicate test, on particle ID array
  // If x is less than node start copy the id to indices_left
  end_pointer = thrust::copy_if(thrust::host,
                                particles->id,
                                particles->id+particles->local_count,
                                particles->x,
                                communication->out_of_bounds.indices_right,
                                GreaterThan(params->node_end_x));
  int oob_right_count = end_pointer - communication->out_of_bounds.indices_right;
  communication->out_of_bounds.particle_count_right = oob_right_count;

  // Remove OOB entries from id array...I wish this was combined with the above
  end_pointer = thrust::remove_if(thrust::host,
                                  particles->id,
                                  particles->id+particles->local_count,
                                  particles->x,
                                  OutsideBounds(params->node_start_x, params->node_end_x));

  // Pack removed particle components
  PackOOBComponents(communication,
                    particles);

  // Set new particle count
  particles->local_count -= (oob_left_count+oob_right_count);

  // Use ID component to reorganize rest of values
  for(int i=0; i<particles->local_count; ++i) {
      const int current_id = particles->id[i];
      // Don't know for sure if testing is worth it...can just copy all of them
      if(current_id != i) {
        CopyParticle(particles, current_id, i);
      }
  }

}
