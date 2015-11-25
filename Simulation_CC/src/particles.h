#pragma once

#include "dimension.h"
#include "vec.h"
#include "hemi/array.h"
#include "hemi/parallel_for.h"
#include "parameters.h"
#include "neighbors.h"
#include "kernels.h"

/**
  @brief Class to handle 2D and 3D SPH particle physics
**/

template<typename Real, Dimension Dim>
class Particles {
public:
  /**
    @brief Constructor: allocates maximum particle storage
  **/
  Particles(const Parameters<Real, Dim>& parameters):
                               m_parameters{parameters},
                               m_neighbors{*this},
                               m_kernels{m_parameters},
                               m_local_count{0},
                               m_max_local_count{m_parameters.GetMaxParticlesLocal()},
                               m_positions{m_max_local_count},
                               m_position_stars{m_max_local_count},
                               m_velocities{m_max_local_count},
                               m_delta_positions{m_max_local_count},
                               m_vorticities{m_max_local_count},
                               m_colors{m_max_local_count},
                               m_densities{m_max_local_count},
                               m_lambdas{m_max_local_count},
                               m_ids{m_max_local_count} {};

  ~Particles()                           = default;
  Particles(const Particles&)            = default;
  Particles& operator=(const Particles&) = default;
  Particles(Particles&&) noexcept        = default;
  Particles& operator=(Particles&&)      = default;

  /**
    @return number of node local particles currently in use
  **/
  std::size_t GetLocalCount() const { return m_local_count; }

  /**
    @return pointer to Vec<> position coordinates array
  **/
  const Vec<Real,Dim>* GetPositionsPointer() const { return m_positions.readOnlyPtr(); }

  /**
    @brief Set number of local particles currently in use
  **/
  void SetLocalCount(size_t new_count) {
    if(new_count < m_max_local_count) {
      m_local_count++;

      m_positions.resize(new_count);
      m_position_stars.resize(new_count);
      m_velocities.resize(new_count);
      m_delta_positions.resize(new_count);
      m_vorticities.resize(new_count);
      m_colors.resize(new_count);
      m_densities.resize(new_count);
      m_lambdas.resize(new_count);
      m_ids.resize(new_count);
    }
    else {
      throw std::runtime_error("too many particles local!");
    }
  }

  /**
    @brief Add fluid particles to current node, filling aabb
  **/
  void AddFluid(const AABB<Real,Dim>& aabb) {
    // |--|--|--|
    // -o--o--o-
    const Vec<Real,Dim> space_counts{Floor((aabb.max - aabb.min) /
                             m_parameters.GetParticleRestSpacing())};

    Vec<int,3> particle_counts(static_cast< Vec<int,Dim> >(space_counts));
    particle_counts.z = particle_counts.z == 0 ? 1 : particle_counts.z;

    std::cout<< particle_counts.x<<","<<particle_counts.y<<","<<particle_counts.z<<std::endl;

    for(int z=0; z < particle_counts[2]; ++z) {
      for(int y=0; y < particle_counts[1]; ++y) {
        for(int x=0; x < particle_counts[0]; ++x) {
          Vec<int,3> int_coord{x,y,z};
          Vec<Real,Dim> coord = static_cast< Vec<Real,Dim> >(int_coord) * m_parameters.GetParticleRestSpacing();
          coord += m_parameters.GetParticleRestSpacing()/static_cast<Real>(2.0);
          const size_t particle_id = this->GetLocalCount();

          m_positions.writeOnlyPtr()[particle_id] = coord;
          m_position_stars.writeOnlyPtr()[particle_id] = coord;
          m_velocities.writeOnlyPtr()[particle_id] = Vec<Real,Dim>{0.0};
          m_ids.writeOnlyPtr()[particle_id] = particle_id;

          this->SetLocalCount(this->GetLocalCount()+1);
        }  // x
      } // y
    } // z
  }

  /**
    @brief Apply external forces to particles
  **/
  void ApplyExternalForces() {
    const Real g = m_parameters.GetGravity();
    const Real dt = m_parameters.GetTimeStep();
    auto velocities = m_velocities.ptr();

    const auto begin = static_cast<std::size_t>(0);
    const auto end = m_local_count;

    hemi::parallel_for(begin, end, [=] HEMI_LAMBDA (std::size_t i) {
      velocities[i].y += g*dt;
    });
  }

  /**
    @brief predict particle positions and apply boundary conditions
  **/
  void PredictPositions() {
    const Real dt = m_parameters.GetTimeStep();
    const auto velocities = m_velocities.readOnlyPtr();
    const auto positions = m_positions.readOnlyPtr();
    auto position_stars = m_position_stars.ptr();
    const auto boundary = m_parameters.GetBoundary();

    const auto begin = static_cast<std::size_t>(0);
    const auto end = m_local_count;

    hemi::parallel_for(begin, end, [=] HEMI_LAMBDA (std::size_t i) {
      position_stars[i] = positions[i] + (velocities[i] * dt);
      ApplyBoundaryConditions(position_stars[i], boundary);
    });
  }

  void ComputeDensities();
  void UpdatePositions();
  void ComputeLambdas();
  void UpdatePositionStars();
  void UpdateDeltaPositions();
  void UpdateVelocities();
  void ApplySurfaceTension();
  void ApplyViscosity();
  void ComputeVorticity();
  void ApplyVorticityConfinement();

private:
  const Parameters<Real,Dim>& m_parameters;
  Neighbors<Real,Dim> m_neighbors;
  Kernels<Real,Dim> m_kernels;
  std::size_t m_local_count;
  const std::size_t m_max_local_count;

  hemi::Array< Vec<Real,Dim> > m_positions;
  hemi::Array< Vec<Real,Dim> > m_position_stars;
  hemi::Array< Vec<Real,Dim> > m_velocities;
  hemi::Array< Vec<Real,Dim> > m_delta_positions;
  hemi::Array< Vec<Real,Dim> > m_vorticities;
  hemi::Array< Vec<Real,Dim> > m_colors;
  hemi::Array< Real > m_densities;
  hemi::Array< Real > m_lambdas;
  hemi::Array< std::size_t > m_ids;
};

/**
  @brief Confine particle position to boundary AABB
**/
template<typename Real, Dimension Dim>
void ApplyBoundaryConditions(Vec<Real,Dim>& position, const AABB<Real,Dim>& boundary) {
  ClampInPlace(position, boundary.min, boundary.max);
}
