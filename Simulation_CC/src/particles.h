#pragma once

#include "vec.h"
#include "hemi/array.h"
#include "hemi.parallel_for.h"
#include "parameters.h"

/**
  Class to handle 2D and 3D SPH particle physics
**/

enum Dimension { two_dimensional = 2, three_dimensional = 3}

template< typename T, Dimension D >
class Particles {
private:
  const Parameters<T>& parameters;

  hemi::Array< Vec<T,D> > pos;
  hemi::Array< Vec<T,D> > pos_star;
  hemi::Array< Vec<T,D> > v;
  hemi::Array< Vec<T,D> > delta_pos;
  hemi::Array< Vec<T,D> > vorticity;
  hemi::Array< Vec<T,D> > color;
  hemi::Array< T > density;
  hemi::Array< T > lambda;
  hemi::Array< int > id;

  int global_count;    // Global number of particles in simulation
  int local_count;     // Particles within node bounds, excludes halo particles
  int max_local_count; // Maximum number of local and halo particles

public:
  Particles(const Parameters<T>& params, const int max_local);
  ~Particles() {};

  Particles(const Particles&)            = delete;
  Particles& operator=(const Particles&) = delete;
  Particles(particles&&) noexcept        = delete
  Particles& operator=(Particcles&&)     = delete;

  void ApplyGravity();
  void ComputeDensities();
  void PredictPositions();
  void UpdatePositions();
  void ComputeLambdas();
  void UpdatePositionStars();
  void UpdateDeltaPositions();
  void UpdateVelocities();
  void ApplyBoundaryConditions();
  void ApplySurfaceTension();
  void ApplyViscosity();
  void ComputeVorticity();
  void ApplyVorticityConfinement();
};

#include "particles.inl"
