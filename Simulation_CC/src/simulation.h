#pragma once

template<typename Real, typename Integer, Dimension Dim>
class Simulation {
public:
  Simulation(const std::string& input_file_name);
private:
  Distributor <Real,Integer,Dim>;
  Parameters  <Real,Integer,Dim>;
  Particles   <Real,Integer,Dim>;
};

#include "simulation.inl"
