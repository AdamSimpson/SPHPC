#pragma once

#include "dimension.h"
#include "vec.h"
#include "array.h"
#include "parameters.h"
#include "neighbors.h"
#include "kernels.h"
#include <thrust/iterator/zip_iterator.h>
#include <thrust/partition.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/for_each.h>

/**
  Class to handle 2D and 3D SPH particle physics
**/

template<typename Real, Dimension Dim>
class Particles {
public:
  /**
    @brief Constructor: allocates maximum particle storage
  **/
  Particles(const Parameters<Real, Dim>& parameters):
                               parameters_{parameters},
                               max_local_count_{parameters.max_particles_local()},
                               neighbors_{parameters},
                               positions_{max_local_count_},
                               position_stars_{max_local_count_},
                               velocities_{max_local_count_},
                               densities_{max_local_count_},
                               lambdas_{max_local_count_},
                               scratch_{max_local_count_} {};

  ~Particles()                           = default;
  Particles(const Particles&)            = default;
  Particles& operator=(const Particles&) = default;
  Particles(Particles&&) noexcept        = default;
  Particles& operator=(Particles&&)      = default;

  /**
    @return number of node/process local particles currently in use
  **/
  std::size_t local_count() const { return positions_.size(); }

  /**
     @return number of particles available
   **/
  std::size_t available() const { return positions_.available(); }

  /**
    @return maximum number of node/process local particles available
  **/
  std::size_t max_local_count() const { return max_local_count_; }

  /**
    @return reference to particle sim::arrays
  **/
  sim::Array<Vec<Real,Dim>>& positions() { return positions_; }
  sim::Array<Vec<Real,Dim>>& position_stars() { return position_stars_; }
  sim::Array<Vec<Real,Dim>>& velocities() { return velocities_; }

  const sim::Array<Vec<Real,Dim>>& positions() const { return positions_; }
  const sim::Array<Vec<Real,Dim>>& position_stars() const { return position_stars_; }
  const sim::Array<Vec<Real,Dim>>& velocities() const { return velocities_; }

  /**
    Remove particles from end of array
  **/
  void remove(std::size_t count) {
    positions_.pop_back(count);
    position_stars_.pop_back(count);
    velocities_.pop_back(count);
    densities_.pop_back(count);
    lambdas_.pop_back(count);
    scratch_.pop_back(count);
  }

  /**
     Add particle to end of array
  **/
  void add(const Vec<Real,Dim>& positions,
           const Vec<Real,Dim>& position_stars,
           const Vec<Real,Dim>& velocities) {

    positions_.push_back(positions);
    position_stars_.push_back(position_stars);
    velocities_.push_back(velocities);

    densities_.push_back(0.0);
    lambdas_.push_back(0.0);
    scratch_.push_back(Vec<Real,Dim>{0.0});
  }

  /**
     Add array of particles to end of array
   **/
  void add(const Vec<Real,Dim>* positions,
           const Vec<Real,Dim>* position_stars,
           const Vec<Real,Dim>* velocities,
           std::size_t count) {

    positions_.push_back(positions, count);
    position_stars_.push_back(position_stars, count);
    velocities_.push_back(velocities, count);

    densities_.push_back(0.0);
    lambdas_.push_back(0.0);
    scratch_.push_back(Vec<Real,Dim>{0.0});
  }

  /**
    @brief Set initial fluid particles to current node, filling 2D aabb
  **/
  std::size_t construct_initial_fluid(const AABB<Real,2>& aabb) {
    // |--|--|--|
    // -o--o--o-
    Vec<std::size_t,Dim> particle_counts = bin_count_in_volume(aabb, parameters_.particle_rest_spacing());

    for(std::size_t y=0; y < particle_counts[1]; ++y) {
      for(std::size_t x=0; x < particle_counts[0]; ++x) {
        Vec<std::size_t,2> int_coord{x,y};
        Vec<Real,Dim> coord = static_cast< Vec<Real,Dim> >(int_coord) *
                               parameters_.particle_rest_spacing()    +
                               aabb.min;
        coord += parameters_.particle_rest_spacing()/static_cast<Real>(2.0);
        auto zero_vector = Vec<Real,Dim>{0.0};
        this->add(coord, coord, zero_vector);
      }  // x
    } // y

    return product(particle_counts);
  }

  /**
    @brief Set initial fluid particles to current node, filling 3D aabb
  **/
  std::size_t construct_initial_fluid(const AABB<Real,3>& aabb) {
    // |--|--|--|
    // -o--o--o-
    Vec<std::size_t,Dim> particle_counts = bin_count_in_volume(aabb, parameters_.particle_rest_spacing());

    std::cout<< particle_counts.x<<","<<particle_counts.y<<","<<particle_counts.z<<std::endl;

    for(std::size_t z=0; z < particle_counts[2]; ++z) {
      for(std::size_t y=0; y < particle_counts[1]; ++y) {
        for(std::size_t x=0; x < particle_counts[0]; ++x) {
          Vec<std::size_t,3> int_coord{x,y,z};
          Vec<Real,Dim> coord = static_cast< Vec<Real,Dim> >(int_coord) *
                                 parameters_.particle_rest_spacing()    +
                                 aabb.min;
          coord += parameters_.particle_rest_spacing()/static_cast<Real>(2.0);
          auto zero_vector = Vec<Real,Dim>{0.0};
          this->add(coord, coord, zero_vector);
        }  // x
      } // y
    } // z
    return product(particle_counts);
  }


  void find_neighbors(IndexSpan to_bin_count, IndexSpan to_fill_count) {
    neighbors_.find(to_bin_count, to_fill_count, position_stars_.data());
  }

  /**
    @brief Apply external forces to particles
  **/
  void apply_external_forces(IndexSpan span) {
    const Real g = parameters_.gravity();
    const Real dt = parameters_.time_step();

    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(begin, end, [=] (std::size_t p) {
      velocities_[p].y += g*dt;
    });
  }

  /**
    @brief predict particle positions and apply boundary conditions
  **/
  void predict_positions(IndexSpan span) {
    const Real dt = parameters_.time_step();

    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(begin, end, [=] (std::size_t p) {
      const auto position_star_p_old = positions_[p];
      auto position_star_p_new = position_star_p_old + (velocities_[p] * dt);
      apply_boundary_conditions(position_star_p_new,
                                parameters_.boundary());
      position_stars_[p] = position_star_p_new;
    });

  }

  void compute_densities(IndexSpan span) {
    const Poly6<Real,Dim> W{parameters_.smoothing_radius()};
    const Real W_0 = W(static_cast<Real>(0.0));

    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(begin, end, [=] (std::size_t p) {
        const Real mass_p = parameters_.rest_mass();
        // Own contribution to density
        Real density = mass_p * W_0;

        // Compiler should hoist this out I ASSume?
        // const auto position_star_p = position_stars_[p];

        for(const std::size_t q : neighbors_[p]) {
          const Real mass_q = parameters_.rest_mass();
          const Real r_mag = magnitude(position_stars_[p] - position_stars_[q]);
          density += mass_q * W(r_mag);
        }

        densities_[p] = density;
    });
  }

  void compute_lambdas(IndexSpan span) {
    const Del_Spikey<Real,Dim> Del_W{parameters_.smoothing_radius()};

    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(begin, end, [=] (std::size_t p) {
      const Real constraint = densities_[p]/parameters_.rest_density() - static_cast<Real>(1.0);
      // Clamp constraint to be positive
      const Real C_p = (constraint < 0.0 ? 0.0 : constraint);
      Real sum_C = 0.0;
      Vec<Real,Dim> sum_gradient{0.0};

      for(const std::size_t q : neighbors_[p]) {
        const Real mass_q = parameters_.rest_mass();
        const Vec<Real,Dim> gradient = mass_q * Del_W(position_stars_[p], position_stars_[q]);
        sum_gradient += gradient;
        // Add k = j contribution
        sum_C += magnitude_squared(gradient);
      }
      // k = i contribution
      sum_C += magnitude_squared(sum_gradient);

      // Having particles of different mass would require 1/m_i term
      // http://mmacklin.com/EG2015PBD.pdf

      sum_C /= (parameters_.rest_density() * parameters_.rest_density());

      lambdas_[p] = -C_p/(sum_C + parameters_.lambda_epsilon());
      });
  }

  void compute_delta_positions(IndexSpan span, const int substep){
    const Del_Spikey<Real,Dim> Del_W{parameters_.smoothing_radius()};

    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(begin, end, [=] (std::size_t p) {
        const Real mass_p = parameters_.rest_mass();

        const Real stiffness = 1.0 - pow(static_cast<Real>(1.0) - parameters_.k_stiff(),
                                         static_cast<Real>(1.0)/static_cast<Real>(substep));

        Vec<Real,Dim> dp{0.0};
        for(const std::size_t q : neighbors_[p]) {
          const Real mass_q = parameters_.rest_mass();
          dp += mass_q * (lambdas_[p] + lambdas_[q]) * Del_W(position_stars_[p],
                                                             position_stars_[q]);
        }

        scratch_[p] = stiffness/parameters_.rest_density() * dp;
      });
  };

  void update_position_stars(IndexSpan span){
    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(begin, end, [=] (std::size_t p) {
        // scratch contains delta positions
        const auto position_star_p_old = position_stars_[p];
        auto position_star_p_new = position_stars_[p] + scratch_[p];
        apply_boundary_conditions(position_star_p_new,
                                  parameters_.boundary());
        position_stars_[p] = position_star_p_new;
      });
  };

  void update_velocities(IndexSpan span){
    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(begin, end, [=] (std::size_t p) {
        Vec<Real,Dim> velocity{(position_stars_[p] - positions_[p]) / parameters_.time_step()};
        clamp_in_place(velocity, static_cast<Real>(-1.0)*parameters_.max_speed(), parameters_.max_speed());
        velocities_[p] = velocity;
    });
  }
 
  void update_positions(IndexSpan span){
    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(begin, end, [=] (std::size_t p) {
        positions_[p] = position_stars_[p];
      });

  }

  void apply_surface_tension(IndexSpan span) {
    const Del_Poly6<Real,Dim> Del_W{parameters_.smoothing_radius()};
    const C_Spline<Real,Dim> C{parameters_.smoothing_radius()};
    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(begin, end, [=] (std::size_t p) {
        Vec<Real,Dim> color{0.0};
        for(const std::size_t q : neighbors_[p]) {
          const Real mass_q = parameters_.rest_mass();
          color += Del_W(position_stars_[p],  position_stars_[q]) * mass_q / densities_[q];
        }
        scratch_[p] = parameters_.smoothing_radius() * color;
      });

    thrust::for_each(begin, end, [=] (std::size_t p) {
        Vec<Real,Dim> surface_tension_force{0.0};
        const Real mass_p = parameters_.rest_mass();
        for(const std::size_t q : neighbors_[p]) {
          const Real mass_q = parameters_.rest_mass();
          const Vec<Real,Dim> r = position_stars_[p] - position_stars_[q];
          const Real r_mag = magnitude(r);
          const Vec<Real,Dim> cohesion_force{-parameters_.gamma() * mass_p * mass_q * C(r_mag) * r_mag};
          const Vec<Real,Dim> curvature_force{-parameters_.gamma() * mass_p * scratch_[p] - scratch_[q]};
          const Real K = 2.0*parameters_.rest_density() / (densities_[p] + densities_[q]);
          surface_tension_force += K * (cohesion_force + curvature_force);
        }

        velocities_[p] += surface_tension_force / densities_[p] * parameters_.time_step();
      });

  }



  void apply_viscosity(IndexSpan span) {
    const Poly6<Real,Dim> W{parameters_.smoothing_radius()};
 
    thrust::counting_iterator<std::size_t> begin(span.begin);
    thrust::counting_iterator<std::size_t> end(span.end);

    thrust::for_each(begin, end, [=] (std::size_t p) {
        Vec<Real,Dim> dv{0.0};
        for(const std::size_t q : neighbors_[p]) {
          const Real mass_q = parameters_.rest_mass();
          dv += (velocities_[q] - velocities_[p]) * W(positions_[p], positions_[q]) * mass_q / densities_[q];
          dv *= parameters_.visc_c();
        }
        velocities_[p] += dv;
      });
  }

private:
  const Parameters<Real,Dim>& parameters_;
  const std::size_t max_local_count_;
  Neighbors<Real,Dim> neighbors_;
 
  sim::Array< Vec<Real,Dim> > positions_;
  sim::Array< Vec<Real,Dim> > position_stars_;
  sim::Array< Vec<Real,Dim> > velocities_;
  sim::Array< Real > densities_;
  sim::Array< Real > lambdas_;

  // Scratch values are used for delta_positions, vorticities, and color field
  // Don't overwrite data you depend on
  sim::Array< Vec<Real,Dim> > scratch_;
};

/**
  @brief Confine new position to boundary AABB
**/
/*
template<typename Real, Dimension Dim>
void apply_boundary_conditions(Vec<Real,Dim>& position,
                               const Vec<Real,Dim>& old_position,
                               const Real time_step,
                               const AABB<Real,Dim>& boundary) {
  // Find effective velocity between old position and new_position
  const Vec<Real,Dim> velocity = (position - old_position) / time_step;

  Real time_to_pen_x;
  if(velocity.x < 0)
    time_to_pen_x = (boundary.min.x - old_position.x)/velocity.x;
  else
    time_to_pen_x = (boundary.max.x - old_position.x)/velocity.x;
  Real time_to_pen_y;
  if(velocity.y < 0)
    time_to_pen_y = (boundary.min.y - old_position.y)/velocity.y;
  else
    time_to_pen_y = (boundary.max.y - old_position.y)/velocity.y;

  if(time_to_pen_x <= time_to_pen_y && (position.x < boundary.min.x || position.x > boundary.max.x)) {
    position = old_position + (velocity * (time_to_pen_x+(float)0.0001));
  }
  else if(time_to_pen_y < time_to_pen_x && (position.y < boundary.min.y || position.y > boundary.max.y)){
    position = old_position + (velocity * (time_to_pen_y+(float)0.0001));
  }

  //  clamp_in_place(position, boundary.min, boundary.max);
}
*/

template<typename Real, Dimension Dim>
void apply_boundary_conditions(Vec<Real,Dim>& position, const AABB<Real,Dim>& boundary) {
  clamp_in_place(position, boundary.min, boundary.max);
}

