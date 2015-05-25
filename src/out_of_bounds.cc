#include "out_of_bounds.h"

#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
extern "C" {
  #define restrict
  #include "simulation.h"
  #include "communication.h"
  #include "particles.h"
  #undef restrict
}
#include <iostream>

struct CheckOOBLeft {
  CheckOOBLeft(double min): x_min(min) {}
  double x_min;

  __host__ __device__
  bool operator()(const double x) {
    return !(x < x_min);
  }
};

struct CheckOOBRight {
  CheckOOBRight(double max): x_max(max) {}
  double x_max;

  __host__ __device__
  bool operator()(const double x) {
    return !(x > x_max);
  }
};

// Remove particles less than the node start and save to packed array;
void PackOOB(struct Params *const params,
             struct Particles *const particles,
             struct Communication *const communication) {

  // Use x array as stencil, for predicate test, on particle ID array
  // If x is less than node start "remove" ID from resulting interator
  // Copy any removed id's to indices_left
  int *end_pointer = thrust::remove_copy_if(particles->id,
                                       particles->id+particles->local_count,
                                       particles->x,
                                       communication->out_of_bounds.indices_left,
                                       CheckOOBLeft(params->node_start_x));
  int oob_left_count = end_pointer - communication->out_of_bounds.indices_left;
  communication->out_of_bounds.particle_count_left = oob_left_count;

  std::cout<<"left_oob: "<<oob_left_count<<std::endl;

  // Use x array as stencil, for predicate test, on particle ID array
  // If x is less than node start "remove" ID from resulting interator
  // Copy any removed id's to indices_left
  end_pointer = thrust::remove_copy_if(particles->id,
                                       particles->id+particles->local_count-oob_left_count,
                                       particles->x,
                                       communication->out_of_bounds.indices_right,
                                       CheckOOBRight(params->node_end_x));
  int oob_right_count = end_pointer - communication->out_of_bounds.indices_right;
  communication->out_of_bounds.particle_count_right = oob_right_count;

  double *packed_send_left = communication->particle_send_buffer;
  double *packed_send_right = communication->particle_send_buffer + oob_left_count*17;

  PackOOBComponents(communication,
                    particles,
                    packed_send_left,
                    packed_send_right);

  // Get new particle count
  particles->local_count -= (oob_left_count+oob_right_count);

  // Use ID component to reorganize rest of values
  for(int i=0; i<particles->local_count; ++i) {
      const int current_id = particles->id[i];
      if(current_id != i)
        CopyParticle(particles, current_id, i);
  }
}
