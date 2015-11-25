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
Vec<Real,Dim> ToRealVec(const std::string& input_string);

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
    this->ReadINI(file_name);
    this->DeriveFromInput();
  };

  ~Parameters()                            = default;
  Parameters(const Parameters&)            = default;
  Parameters& operator=(const Parameters&) = default;
  Parameters(Parameters&&) noexcept        = default;
  Parameters& operator=(Parameters&&)      = default;

  /**
    @brief Read file_name .INI containing parameters
  **/
  void ReadINI(const std::string& file_name) {
    boost::property_tree::ptree property_tree;
    boost::property_tree::ini_parser::read_ini(file_name, property_tree);

    m_time_step_count = property_tree.get<std::size_t>("SimParameters.number_steps");
    m_time_step = property_tree.get<Real>("SimParameters.time_step");
    m_global_particle_count = property_tree.get<std::size_t>("SimParameters.global_particle_count");
    m_max_particles_local = property_tree.get<std::size_t>("SimParameters.max_particles_local");
    m_g = property_tree.get<Real>("PhysicalParameters.g");
    m_c = property_tree.get<Real>("PhysicalParameters.c");
    m_k = property_tree.get<Real>("PhysicalParameters.k");

    m_boundary.min = ToRealVec<Real,Dim>(property_tree.get<std::string>("Boundary.min"));
    m_boundary.max = ToRealVec<Real,Dim>(property_tree.get<std::string>("Boundary.max"));

    m_initial_fluid.min = ToRealVec<Real,Dim>(property_tree.get<std::string>("InitialFluid.min"));
    m_initial_fluid.max = ToRealVec<Real,Dim>(property_tree.get<std::string>("InitialFluid.max"));
  }

  /**
    @brief Derive additional parameters from .INI parameters required for simulation
  **/
  void DeriveFromInput() {
    m_particle_rest_spacing = pow(m_initial_fluid.Volume() / m_global_particle_count, 1.0/Dim);
    m_particle_radius = m_particle_rest_spacing/2.0;
    m_smoothing_radius = 2.0*m_particle_rest_spacing;
  }

  /**
    @return maximum node level particle count
  **/
  std::size_t GetMaxParticlesLocal() const {
    return m_max_particles_local;
  }

  /**
    @return maximum global level particle count
  **/
  std::size_t GetGlobalParticleCount() const {
    return m_global_particle_count;
  }

  /**
    Set maximum global level particle count
  **/
  void SetGlobalParticleCount(std::size_t global_count) {
    m_global_particle_count = global_count;
  }

  /**
    @return Initial fluid configuration as AABB
  **/
  AABB<Real,Dim> GetInitialFluid() const {
    return m_initial_fluid;
  }

  /**
    @return Global boundary as AABB
  **/
  AABB<Real,Dim> GetBoundary() const {
    return m_boundary;
  }

  /**
    @return Rest spacing of particles
  **/
  Real GetParticleRestSpacing() const {
    return m_particle_rest_spacing;
  }

  /**
    @return Gravitational acceleration magnitude
  **/
  Real GetGravity() const {
    return m_g;
  }

  /**
    @return number of simulation steps to take
  **/
  std::size_t GetTimeStepCount() const {
    return m_time_step_count;
  }

  /**
    @return time step size
  **/
  Real GetTimeStep() const {
    return m_time_step;
  }

private:
  std::size_t m_max_particles_local;
  std::size_t m_global_particle_count;
  std::size_t m_time_step_count;
  Real m_particle_rest_spacing;
  Real m_particle_radius;
  Real m_smoothing_radius;
  Real m_g;
  Real m_c;
  Real m_k;
  Real m_time_step;
  AABB<Real,Dim> m_boundary;
  AABB<Real,Dim> m_initial_fluid;
};

/**
  @brief Returns a Vec<> from comma seperated input string
**/
template<typename Real, Dimension Dim>
Vec<Real,Dim> ToRealVec(const std::string& input_string) {
  Vec<Real,Dim> result;
  std::stringstream ss(input_string);
  std::string item;

  for(int i=0; i<Dim; i++) {
    std::getline(ss, item, ',');
    boost::algorithm::trim(item);
    result[i] = (boost::lexical_cast<Real>(item));
  }

  return result;
}
