#include "input_parser.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
extern "C"
{
#include "simulation.h"
}

void set_parameters(Params *const parameters) {
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini("sim-config.ini", pt);
  std::cout << pt.get<std::string>("SimParameters.number_steps") << std::endl;
  std::cout << pt.get<std::string>("SimParameters.time_step") << std::endl;
  std::cout << pt.get<std::string>("SimParameters.number_particles") << std::endl;
}
