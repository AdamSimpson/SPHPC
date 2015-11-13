#pragma once

#include <string>
#include "dimension.h"
#include "distributor.h"
#include "parameters.h"
#include "particles.h"

template<typename Real, Dimension Dim>
class Simulation {
public:
  /**
    Simulation constructor
    If constructor fails distributor is likely to abort application
    instead of throwing exception as MPI is finalized with exception
  **/
  Simulation(const std::string& input_file_name) :
    parameters{input_file_name},
    particles{parameters} {};

    ~Simulation()                            = default;
    Simulation(const Simulation&)            = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) noexcept        = delete;
    Simulation& operator=(Simulation&&)      = delete;

private:
  Distributor <Real,Dim> distributor;
  Parameters  <Real,Dim> parameters;
  Particles   <Real,Dim> particles;
};
