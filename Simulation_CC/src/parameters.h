#pragma once

#include "dimension.h"
#include "vec.h"
#include "aabb.h"
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

template<typename Real, Dimension Dim>
Vec<Real,Dim> ToRealVec(const std::string& input_string);

template<typename Real, Dimension Dim>
class Parameters {
public:
  /**
  Parameters constructor:
    Handle exceptions as they will fail to propogate to main()
  **/
  Parameters(const std::string& file_name) {
    ReadParameters(file_name);
  };

  ~Parameters()                            = default;
  Parameters(const Parameters&)            = default;
  Parameters& operator=(const Parameters&) = default;
  Parameters(Parameters&&) noexcept        = default;
  Parameters& operator=(Parameters&&)      = default;

  /**
    Read INI file containing parameter files
  **/
  void ReadParameters(const std::string& file_name) {
    boost::property_tree::ptree property_tree;
    boost::property_tree::ini_parser::read_ini(file_name, property_tree);

    time_step_count = property_tree.get<std::size_t>("SimParameters.number_steps");
    time_step = property_tree.get<Real>("SimParameters.time_step");
    global_particle_count = property_tree.get<std::size_t>("SimParameters.global_particles_local");
    max_particles_local = property_tree.get<std::size_t>("SimParameters.max_particles_local");
    g = property_tree.get<Real>("PhysicalParameters.g");
    c = property_tree.get<Real>("PhysicalParameters.c");
    k = property_tree.get<Real>("PhysicalParameters.k");

    boundary.min = ToRealVec<Real,Dim>(property_tree.get<std::string>("Boundary.min"));
    boundary.max = ToRealVec<Real,Dim>(property_tree.get<std::string>("Boundary.max"));

    initial_fluid.min = ToRealVec<Real,Dim>(property_tree.get<std::string>("InitialFluid.min"));
    initial_fluid.max = ToRealVec<Real,Dim>(property_tree.get<std::string>("InitialFluid.min"));
  }

  std::size_t GetMaxParticlesLocal() const {
    return max_particles_local;
  };

private:
  std::size_t max_particles_local;
  std::size_t global_particle_count;
  std::size_t time_step_count;
  Real g;
  Real c;
  Real k;
  Real time_step;
  AABB<Real,Dim> boundary;
  AABB<Real,Dim> initial_fluid;
};

// Returns a Vec<> from comma seperated input string
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
