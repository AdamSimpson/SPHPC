#include "particles.hpp"
#include "hemi/hemi.h"
#include <iostream>

template<typename T, int dim>
Particles::Particles(const int max_local_) : star{hemi::Array< Vector<T, dim> >(max_local_)},
                                             pos{hemi::Array< Vector<T, dim> >(max_local_)},
                                             v{hemi::Array< Vector<T, dim> >(max_local_)},
                                             dp{hemi::Array< Vector<T, dim> >(max_local_)},
                                             w{hemi::Array< Vector<T, dim> >(max_local_)},
                                             color{hemi::Array< Vector<T, dim> >(max_local_)},
                                             density{hemi::Array<T>(max_local_)},
                                             density{hemi::Array<T>(max_local_)},
                                             density{hemi::Array<int>(max_local_)} {

  ConstructFluidVolume(particles, params, fluid_volume_initial);

  rest_density = 1000.0;
  rest_mass = params->fluid_volume * rest_density / (double)global_count;
  // For some reason the initial density estimate is slightly too high
  // Perhaps fluid_volume of constructed fluid is incorrect?
  std::cout<<"Fudging rest mass\n";
  rest_mass *= .991;

  // Initialize particle values`
  for (int i=0; i<local_count; ++i) {
    v_x[i] = 0.0;
    v_y[i] = 0.0;
    v_z[i] = 0.0;
    x_star[i] = particles->x[i];
    y_star[i] = particles->y[i];
    z_star[i] = particles->z[i];
    density[i] = particles->rest_density;
  }

  halo_count_left = 0;
  halo_count_right = 0;
}
