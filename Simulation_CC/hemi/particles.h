#ifndef SPH_SRC_PARTICLES_H_
#define SPH_SRC_PARTICLES_H_

#include "vec.h"
#include "hemi/array.h"
#include "hemi/parallel_for.h"
#include "parameters.h"

/**
Class to handle 2D and 3D SPH particles
**/

enum Dimension { two_dimensional=2, three_dimensional=3 };

template<typename TReal, Dimension TDim>
class Particles {
private:
  const Parameters<TReal>& parameters;
  hemi::Array< Vector<TReal, TDim> > pos;
  hemi::Array< Vector<TReal, TDim> > v;
  hemi::Array<TReal> density;

  int global_count; // Global number of particles in simulation
  int local_count; // Particles within node bounds, excludes halo particles
  int max_count_local;

  // Implimentation of template functions
  #include "particles.hxx"

};

/**
Explicitly instantiate valid particle templates
**/
template class Particles<float,  two_dimensional>;
template class Particles<float,  three_dimensional>;
template class Particles<double, two_dimensional>;
template class Particles<double, three_dimensional>;
#endif
