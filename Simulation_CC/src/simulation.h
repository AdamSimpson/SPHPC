#pragma once

#include <string>
#include "dimension.h"
#include "distributor.h"
#include "parameters.h"
#include "particles.h"
#include "aabb.h"

template<typename Real, Dimension Dim>
class Simulation {
public:
  /**
    Simulation constructor
    Reads in parameters and constructs initial fluid
  **/
  Simulation(const std::string& input_file_name) :
    distributor{parameters, particles},
    parameters{input_file_name},
    particles{parameters} {
      distributor.Initilize();
      this->AddFluid(parameters.GetInitialFluid());
    };

    ~Simulation()                            = default;
    Simulation(const Simulation&)            = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) noexcept        = delete;
    Simulation& operator=(Simulation&&)      = delete;

    void AddFluid(const AABB<Real,Dim>& aabb) {
      distributor.DistributeFluid(parameters.GetInitialFluid());
    }

private:
  Distributor <Real,Dim> distributor;
  Parameters  <Real,Dim> parameters;
  Particles   <Real,Dim> particles;
};
