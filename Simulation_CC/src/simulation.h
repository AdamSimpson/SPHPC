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
    distributor_{parameters_, particles_},
    parameters_{input_file_name},
    particles_{parameters_},
    writer_{"adios_writer.xml", distributor_, particles_} {
      distributor_.initilize();
      this->construct_initial_fluid(parameters_.initial_fluid());
    }

    ~Simulation()                            = default;
    Simulation(const Simulation&)            = delete;
    Simulation& operator=(const Simulation&) = delete;
    Simulation(Simulation&&) noexcept        = delete;
    Simulation& operator=(Simulation&&)      = delete;

    void construct_initial_fluid(const AABB<Real,Dim>& aabb) {
      distributor_.distribute_initial_fluid(aabb);
    }

    void write_particles() {
      writer_.write_particles();
    }

    void print_global_particle_count() {
      std::cout<<distributor_.global_count()<<std::endl;
    }

    std::size_t time_step_count() {
      return parameters_.time_step_count();
    }

    void advance_particles() {
      std::cout<<"advancing "<<std::endl;
      particles_.apply_external_forces(distributor_.resident_span());
      particles_.predict_positions(distributor_.resident_span());

      // @todo Balance nodes

      //      distributor_.domain_sync();

      particles_.find_neighbors(distributor_.local_span(),
                                distributor_.resident_span());

      for(int sub=0; sub<parameters_.solve_step_count(); sub++) {
        particles_.compute_densities(distributor_.local_span());
        // update halo densities
        particles_.compute_lambdas(distributor_.local_span());
        // update halo lambdas
        particles_.compute_delta_positions(distributor_.local_span(), sub);
        // update halo deltas
        particles_.update_position_stars(distributor_.local_span());
        // update position star halos?
      }

      particles_.update_velocities(distributor_.local_span());

      particles_.apply_surface_tension(distributor_.local_span());

      particles_.apply_viscosity(distributor_.local_span());

      // viscosity
      // vorticity
      // vorticity confinement

      particles_.update_positions(distributor_.local_span());
    }


private:
  Distributor <Real,Dim> distributor_;
  Parameters  <Real,Dim> parameters_;
  Particles   <Real,Dim> particles_;
  AdiosWriter <Real,Dim> writer_;
};
