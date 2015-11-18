#pragma once

#include "dimension.h"
#include "vec.h"
#include "hemi/array.h"
#include "hemi/parallel_for.h"
#include "parameters.h"
#include "neighbors.h"
#include "kernels.h"

/**
  Class to handle 2D and 3D SPH particle physics
**/

template<typename Real, Dimension Dim>
class Particles {
public:
  Particles(const Parameters<Real, Dim>& params):
                               parameters{params},
                               neighbors{*this},
                               kernels{params},
                               max_local_count{parameters.GetMaxParticlesLocal()},
                               pos{max_local_count},
                               pos_star{max_local_count},
                               v{max_local_count},
                               delta_pos{max_local_count},
                               vorticity{max_local_count},
                               color{max_local_count},
                               density{max_local_count},
                               lambda{max_local_count},
                               id{max_local_count} {};

  ~Particles()                           = default;
  Particles(const Particles&)            = default;
  Particles& operator=(const Particles&) = default;
  Particles(Particles&&) noexcept        = default;
  Particles& operator=(Particles&&)      = default;

  std::size_t GetLocalCount() const { return local_count; }
  /**
    Set number of local particles
  **/
  void SetLocalCount(size_t new_count) {
    if(new_count < max_local_count) {
      local_count++;

      pos.resize(local_count);
      pos_star.resize(local_count);
      v.resize(local_count);
      delta_pos.resize(local_count);
      vorticity.resize(local_count);
      color.resize(local_count);
      density.resize(local_count);
      lambda.resize(local_count);
      id.resize(local_count);
    }
    else {
      throw std::string("Too many particles!\n");
    }
  }

  /**
    Add fluid particles to current node
  **/
  void AddFluid(const AABB<Real,Dim>& aabb) {
    // o--o--o--o : space_counts == 3, particle_counts = 4
    const Vec<Real,Dim> space_counts{Floor((aabb.max - aabb.min) /
                             parameters.GetParticleSpacing())};

    Vec<int,3> particle_counts(static_cast< Vec<int,Dim> >(space_counts));
    particle_counts += 1;

    for(int z=0; z < particle_counts[2]; ++z) {
      for(int y=0; y < particle_counts[1]; ++y) {
        for(int x=0; x < particle_counts[0]; ++x) {
          const Vec<int,3> int_coord{x+1,y+1,z+1};
          const Vec<Real,Dim> coord = static_cast< Vec<Real,Dim> >(int_coord) * parameters.GetParticleSpacing();
          const size_t particle_id = this->GetLocalCount();
/*
          pos[particle_id] = coord;
          pos_star[particle_id] = coord;
          v[particle_id] = Vec<Real,Dim>{0.0};
          id[particle_id] = particle_id;
*/
          this->SetLocalCount(this->GetLocalCount()+1);
        }  // x
      } // y
    } // z
  }

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
  Kernels<Real,Dim> kernels;
  std::size_t local_count;
  const std::size_t max_local_count;

  hemi::Array< Vec<Real,Dim> > pos;
  hemi::Array< Vec<Real,Dim> > pos_star;
  hemi::Array< Vec<Real,Dim> > v;
  hemi::Array< Vec<Real,Dim> > delta_pos;
  hemi::Array< Vec<Real,Dim> > vorticity;
  hemi::Array< Vec<Real,Dim> > color;
  hemi::Array< Real > density;
  hemi::Array< Real > lambda;
  hemi::Array< std::size_t > id;
};
