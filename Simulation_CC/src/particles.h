#pragma once

#include "vec.h"
#include "hemi/array.h"
#include "hemi.parallel_for.h"
#include "parameters.h"

/**
  Class to handle 2D and 3D SPH particle physics
**/

enum Dimension { two_dimensional = 2, three_dimensional = 3}

template<typename Real, typename Integer, Dimension Dim>
class Particles {
private:
  const Parameters<Real,Integer,Dimension>& parameters;

  typedef Vector_t Vec<Real,Integer,Dimension>
  hemi::Array< Vector_t > pos;
  hemi::Array< Vector_t > pos_star;
  hemi::Array< Vector_t > v;
  hemi::Array< Vector_t > delta_pos;
  hemi::Array< Vector_t > vorticity;
  hemi::Array< Vector_t > color;
  hemi::Array< Real > density;
  hemi::Array< Real > lambda;
  hemi::Array< Integer > id;

  Integer global_count;    // Global number of particles in simulation
  Integer local_count;     // Particles within node bounds, excludes halo particles
  Integer max_local_count; // Maximum number of local and halo particles

public:
  Particles(const Parameters<Real>& params, const int max_local);
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
