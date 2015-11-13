#pragma once

#include "dimension.h"
#include "vec.h"
#include "hemi/array.h"
#include "hemi/parallel_for.h"
#include "parameters.h"
#include "neighbors.h"

/**
  Class to handle 2D and 3D SPH particle physics
**/

template<typename Real, Dimension Dim>
class Particles {
public:
  Particles(const Parameters<Real, Dim>& params):
                               parameters{params},
                               neighbors{*this},
                               pos{parameters.GetMaxParticlesLocal()},
                               pos_star{parameters.GetMaxParticlesLocal()},
                               v{parameters.GetMaxParticlesLocal()},
                               delta_pos{parameters.GetMaxParticlesLocal()},
                               vorticity{parameters.GetMaxParticlesLocal()},
                               color{parameters.GetMaxParticlesLocal()},
                               density{parameters.GetMaxParticlesLocal()},
                               lambda{parameters.GetMaxParticlesLocal()},
                               id{parameters.GetMaxParticlesLocal()} {};

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
  const Parameters<Real,Dim>& parameters;
  Neighbors<Real,Dim> neighbors;

  hemi::Array< Vec<Real,Dim> > pos;
  hemi::Array< Vec<Real,Dim> > pos_star;
  hemi::Array< Vec<Real,Dim> > v;
  hemi::Array< Vec<Real,Dim> > delta_pos;
  hemi::Array< Vec<Real,Dim> > vorticity;
  hemi::Array< Vec<Real,Dim> > color;
  hemi::Array< Real > density;
  hemi::Array< Real > lambda;
  hemi::Array< std::size_t > id;

  std::size_t global_count;    // Global number of particles in simulation
  std::size_t local_count;     // Particles within node bounds, excludes halo particles
  std::size_t max_local_count; // Maximum number of local and halo particles
};
