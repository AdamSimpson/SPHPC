#include "parameters.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <iostream>

void Parameters::ReadParameters() {
  boost::property_tree::ptree property_tree;

  try {
    boost::property_tree::ini_parser::read_ini(this->ini_file_name, property_tree);

    // Camera Parameters
    this->camera.view_up     = ToDouble3(property_tree.get<std::string>("Camera.view_up"));
    this->camera.position    = ToDouble3(property_tree.get<std::string>("Camera.position"));
    this->camera.focal_point = ToDouble3(property_tree.get<std::string>("Camera.focal_point"));

    // Boundary Parameters
    this->boundary.min_coord = ToDouble3(property_tree.get<std::string>("Boundary.min_coord"));
    this->boundary.max_coord = ToDouble3(property_tree.get<std::string>("Boundary.max_coord"));

    // Obstacle Parameters
    this->obstacle.OBJ_file_name = property_tree.get<std::string>("Obstacle.OBJ_file_name");
    this->obstacle.MTL_file_name = property_tree.get<std::string>("Obstacle.MTL_file_name");
    this->obstacle.min_coord     = ToDouble3(property_tree.get<std::string>("Obstacle.min_coord"));
    this->obstacle.max_x         = property_tree.get<double>("Obstacle.max_x");

    this->input_file_path = property_tree.get<std::string>("Input.file_path");

  } catch(std::exception const& exception) {
      std::cout << "Aborting: " << exception.what() << std::endl;
      exit(-1);
  }
}

// Returns a double3 from comma seperated input string
double3 ToDouble3(const std::string input_string) {
  double3 result;
  std::stringstream ss(input_string);
  std::string item;
  std::getline(ss, item, ',');
  boost::algorithm::trim(item);
  result.x = (boost::lexical_cast<double>(item));
  std::getline(ss, item, ',');
  boost::algorithm::trim(item);
  result.y = (boost::lexical_cast<double>(item));
  std::getline(ss, item);
  boost::algorithm::trim(item);
  result.z = (boost::lexical_cast<double>(item));

  return result;
}
