#pragma once

#include "dimension.h"
#include "vec.h"
#include "aabb.h"
#include <cmath>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

// Forward declaration
template<typename Real, Dimension Dim>
Vec<Real,Dim> to_real_vec(const std::string& input_string);

/**
  @brief Simulation wide tunable parameters
**/
template<typename Real, Dimension Dim>
class Parameters {
public:

  /**
    @brief  Construct initial parameters from file_name .INI file
  **/
  Parameters(const std::string& file_name) {
    this->read_INI(file_name);
    this->derive_from_input();
  };

  ~Parameters()                            = default;
  Parameters(const Parameters&)            = default;
  Parameters& operator=(const Parameters&) = default;
  Parameters(Parameters&&) noexcept        = default;
  Parameters& operator=(Parameters&&)      = default;

  /**
    @brief Read file_name .INI containing parameters
  **/
  void read_INI(const std::string& file_name) {
    boost::property_tree::ptree property_tree;
    boost::property_tree::ini_parser::read_ini(file_name, property_tree);

    time_step_count_ = property_tree.get<std::size_t>("SimParameters.number_steps");
    solve_step_count_ = property_tree.get<std::size_t>("SimParameters.number_solve_steps");
    time_step_ = property_tree.get<Real>("SimParameters.time_step");
    initial_global_particle_count_ = property_tree.get<std::size_t>("SimParameters.global_particle_count");
    max_particles_local_ = property_tree.get<std::size_t>("SimParameters.max_particles_local");
    g_ = property_tree.get<Real>("PhysicalParameters.g");
    c_ = property_tree.get<Real>("PhysicalParameters.c");
    gamma_ = property_tree.get<Real>("PhysicalParameters.gamma");
    visc_c_ = property_tree.get<Real>("PhysicalParameters.visc_c");
    lambda_epsilon_ = property_tree.get<Real>("PhysicalParameters.lambda_epsilon");
    k_stiff_ = property_tree.get<Real>("PhysicalParameters.k_stiff");
    rest_density_ = property_tree.get<Real>("PhysicalParameters.density");

    boundary_.min = to_real_vec<Real,Dim>(property_tree.get<std::string>("Boundary.min"));
    boundary_.max = to_real_vec<Real,Dim>(property_tree.get<std::string>("Boundary.max"));

    initial_fluid_.min = to_real_vec<Real,Dim>(property_tree.get<std::string>("InitialFluid.min"));
    initial_fluid_.max = to_real_vec<Real,Dim>(property_tree.get<std::string>("InitialFluid.max"));
  }

  /**
    @brief Derive additional parameters from .INI parameters required for simulation
  **/
  void derive_from_input() {
    particle_rest_spacing_ = pow(initial_fluid_.volume() / initial_global_particle_count_, 1.0/Dim);
    particle_radius_ = particle_rest_spacing_/2.0;
    smoothing_radius_ = 2.0*particle_rest_spacing_;

    Vec<std::size_t,Dim> particle_counts = bin_count_in_volume(initial_fluid_, particle_rest_spacing_);
    std::size_t particle_count = product(particle_counts);
    rest_mass_ = initial_fluid_.volume() * rest_density_ / particle_count;

    max_speed_ = smoothing_radius_ / time_step_;
  }

  /**
    @return maximum node level particle count
  **/
  std::size_t max_particles_local() const {
    return max_particles_local_;
  }

  /**
    Set Initial requested global particle count
  **/
  void set_initial_global_particle_count(std::size_t global_count) {
    initial_global_particle_count_ = global_count;
  }

  /**
    @return Initial fluid configuration as AABB
  **/
  AABB<Real,Dim> initial_fluid() const {
    return initial_fluid_;
  }

  /**
    @return Global boundary as AABB
  **/
  AABB<Real,Dim> boundary() const {
    return boundary_;
  }

  /**
    @return Rest spacing of particles
  **/
  Real particle_rest_spacing() const {
    return particle_rest_spacing_;
  }

  /**
    @return smoothing radius for particles
  **/
  Real smoothing_radius() const {
    return smoothing_radius_;
  }

  /**
    @return Gravitational acceleration magnitude
  **/
  Real gravity() const {
    return g_;
  }

  /**
    @return number of simulation steps to take
  **/
  std::size_t time_step_count() const {
    return time_step_count_;
  }

  /**
    @return time step size
  **/
  Real time_step() const {
    return time_step_;
  }

  /**
     @return solve step count
   **/
  std::size_t solve_step_count() const {
    return solve_step_count_;
  }

  Real rest_mass() const {
    return rest_mass_;
  }

  Real rest_density() const {
    return rest_density_;
  }

  Real lambda_epsilon() const {
    return lambda_epsilon_;
  }

  Real k_stiff() const {
    return k_stiff_;
  }

  Real max_speed() const {
    return max_speed_;
  }

  Real gamma() const {
    return gamma_;
  }

  Real visc_c() const {
    return visc_c_;
  }

private:
  std::size_t max_particles_local_;
  std::size_t initial_global_particle_count_;
  std::size_t time_step_count_;
  std::size_t solve_step_count_;
  Real particle_rest_spacing_;
  Real particle_radius_;
  Real smoothing_radius_;
  Real rest_density_;
  Real rest_mass_;
  Real g_;
  Real c_;
  Real gamma_;
  Real lambda_epsilon_;
  Real k_stiff_;
  Real visc_c_;
  Real time_step_;
  Real max_speed_;
  AABB<Real,Dim> boundary_;
  AABB<Real,Dim> initial_fluid_;
};

/**
  @brief Returns a Vec<> from comma seperated input string
**/
template<typename Real, Dimension Dim>
Vec<Real,Dim> to_real_vec(const std::string& input_string) {
  Vec<Real,Dim> result;
  std::stringstream ss(input_string);
  std::string item;

  for(int i=0; i<Dim; ++i) {
    std::getline(ss, item, ',');
    boost::algorithm::trim(item);
    result[i] = (boost::lexical_cast<Real>(item));
  }

  return result;
}
