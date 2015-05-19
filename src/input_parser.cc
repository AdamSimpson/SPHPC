#include "input_parser.h"

extern "C" {
#include "simulation.h"
#include "geometry.h"
}

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iostream>

void ReadParameters(struct Params *const parameters,
                    struct AABB *boundary_volume,
                    struct AABB *initial_fluid_volume) {
  boost::property_tree::ptree property_tree;

  try {
    boost::property_tree::ini_parser::read_ini(
        "./sim-config.ini", property_tree);

    parameters->number_steps = property_tree.get<int>(
        "SimParameters.number_steps");
    parameters->time_step = property_tree.get<double>(
        "SimParameters.time_step");
    parameters->number_particles_global = property_tree.get<int>(
        "SimParameters.number_particles");
    parameters->g = property_tree.get<double>("PhysicalParameters.g");
    parameters->c = property_tree.get<double>("PhysicalParameters.c");
    parameters->k = property_tree.get<double>("PhysicalParameters.k");

    boundary_volume->min_x = property_tree.get<double>("BoundaryVolume.min_x");
    boundary_volume->min_y = property_tree.get<double>("BoundaryVolume.min_y");
    boundary_volume->min_z = property_tree.get<double>("BoundaryVolume.min_z");
    boundary_volume->max_x = property_tree.get<double>("BoundaryVolume.max_x");
    boundary_volume->max_y = property_tree.get<double>("BoundaryVolume.max_y");
    boundary_volume->max_z = property_tree.get<double>("BoundaryVolume.max_z");

    initial_fluid_volume->min_x = property_tree.get<double>(
        "InitialFluidVolume.min_x");
    initial_fluid_volume->min_y = property_tree.get<double>(
        "InitialFluidVolume.min_y");
    initial_fluid_volume->min_z = property_tree.get<double>(
        "InitialFluidVolume.min_z");
    initial_fluid_volume->max_x = property_tree.get<double>(
        "InitialFluidVolume.max_x");
    initial_fluid_volume->max_y = property_tree.get<double>(
        "InitialFluidVolume.max_y");
    initial_fluid_volume->max_z = property_tree.get<double>(
        "InitialFluidVolume.max_z");
  } catch(std::exception const& exception) {
      std::cout << "Aborting: " << exception.what() << std::endl;
      exit(-1);
  }

}
