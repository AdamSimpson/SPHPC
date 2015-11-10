#pragma once

#include <string>
#include "dimension.h"
#include "distributor.h"
#include "parameters.h"
#include "particles.h"

template<typename Real, typename Integer, Dimension Dim>
class Simulation {
public:
  Simulation(const std::string& input_file_name):
    distributor{Distributor<Real,Integer,Dim>()},
    parameters{Parameters<Real,Integer,Dim>(input_file_name)},
    particles{Particles<Real,Integer,Dim>(parameters)} {/* Set distributor parameters*/};

    ~Simulation()                            =default;
    Simulation(const Simulation&)            = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) noexcept        = delete;
    Simulation& operator=(Simulation&&)      = delete;

private:
  Distributor <Real,Integer,Dim> distributor;
  Parameters  <Real,Integer,Dim> parameters;
  Particles   <Real,Integer,Dim> particles;
};
