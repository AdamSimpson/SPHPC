#pragma once

#include <string>
#include "dimension.h"
#include "distributor.h"
#include "parameters.h"
#include "particles.h"

template<typename Real, typename Integer, Dimension Dim>
class Simulation {
public:
  Simulation(const std::string& input_file_name);
private:
  Distributor <Real,Integer,Dim> distributor;
  Parameters  <Real,Integer,Dim> parameters;
  Particles   <Real,Integer,Dim> particles;
};

#include "simulation.inl"
