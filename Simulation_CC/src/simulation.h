#pragma once

#include <string>
#include "dimension.h"
#include "distributor.h"
#include "parameters.h"
#include "particles.h"
#include "adios_writer.h"
#include "aabb.h"

template<typename Real, Dimension Dim>
class Simulation {
public:
  /**
    Simulation constructor
    Reads in parameters and constructs initial fluid
    responsible for coordinating simulation components
  **/
  Simulation(const std::string& input_file_name) :
    m_distributor{m_parameters, m_particles},
    m_parameters{input_file_name},
    m_particles{m_parameters},
    m_writer{"adios_writer.xml", m_distributor, m_particles, m_parameters} {
      m_distributor.Initilize();
      this->AddFluid(m_parameters.GetInitialFluid());
    }

    ~Simulation()                            = default;
    Simulation(const Simulation&)            = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) noexcept        = delete;
    Simulation& operator=(Simulation&&)      = delete;

    void AddFluid(const AABB<Real,Dim>& aabb) {
      m_distributor.DistributeFluid(m_parameters.GetInitialFluid());
      m_parameters.SetGlobalParticleCount(m_distributor.GetGlobalParticleCount());
    }

    void WriteParticles() {
      m_writer.WriteParticles();
    }

    void PrintGlobalParticleCount() {
      std::cout<<m_distributor.GetGlobalParticleCount()<<std::endl;
    }

    std::size_t GetTimeStepCount() {
      return m_parameters.GetTimeStepCount();
    }

    void AdvanceParticles() {
      m_particles.ApplyExternalForces();
      m_particles.PredictPositions();

      m_distributor.ExchangeOOB();
      m_distributor.ExchangeHalo();
    }

private:
  Distributor <Real,Dim> m_distributor;
  Parameters  <Real,Dim> m_parameters;
  Particles   <Real,Dim> m_particles;
  AdiosWriter <Real,Dim> m_writer;
};
