#pragma once

#include "dimension.h"
#include "vec.h"
#include "hemi/array.h"
#include "hemi/parallel_for.h"
#include "parameters.h"

/**
  Class to handle 2D and 3D SPH particle physics
**/

template<typename Real, typename Integer, Dimension Dim>
class Particles {
public:
  Particles(const Parameters<Real, Integer, Dim>& params);
  ~Particles()                           = default;
  Particles(const Particles&)            = default;
  Particles& operator=(const Particles&) = default;
  Particles(Particles&&) noexcept        = default;
  Particles& operator=(Particles&&)      = default;

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

private:
  const Parameters<Real,Integer,Dim>& parameters;

  hemi::Array< Vec<Real,Dim> > pos;
  hemi::Array< Vec<Real,Dim> > pos_star;
  hemi::Array< Vec<Real,Dim> > v;
  hemi::Array< Vec<Real,Dim> > delta_pos;
  hemi::Array< Vec<Real,Dim> > vorticity;
  hemi::Array< Vec<Real,Dim> > color;
  hemi::Array< Real > density;
  hemi::Array< Real > lambda;
  hemi::Array< Integer > id;

  Integer global_count;    // Global number of particles in simulation
  Integer local_count;     // Particles within node bounds, excludes halo particles
  Integer max_local_count; // Maximum number of local and halo particles
};

#include "particles.inl"
