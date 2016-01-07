#ifndef SPH_SRC_PARTICLES_H_
#define SPH_SRC_PARTICLES_H_

#include "vector.hpp"
#include "hemi/hemi.h"

/**
Class to handle SPH particles
**/

template<typename T, int dim>
class Particles {
  hemi::Array< Vector<T, dim> > star;
  hemi::Array< Vector<T, dim> > pos;
  hemi::Array< Vector<T, dim> > v;
  hemi::Array< Vector<T, dim> > dp;
  hemi::Array< Vector<T, dim> > w;
  hemi::Array< Vector<T, dim> > color;
  hemi::Array<T> density;
  hemi::Array<T> lambda;
  hemi::Array<int> id;

  T rest_density;
  T rest_mass;
  T rest_radius;
  int global_count; // Global number of particles in simulation
  int local_count; // Particles within node bounds, excludes halo particles
  int max_local;  // Maximum number of local and halo particles
  int halo_count_left;
  int halo_count_right;

  Particles(const int max_local_);
};

/**
Explicitly instantiate template
**/
template class Particles<float,  2>;
template class Particles<float,  3>;
template class Particles<double, 2>;
template class Particles<double, 3>;

#endif
