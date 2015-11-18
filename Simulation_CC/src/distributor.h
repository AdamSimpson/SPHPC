#pragma once

#include <boost/mpi.hpp>
#include <utility>
#include <string>
#include "utility_math.h"
#include "particles.h"
#include "parameters.h"

template<typename Real, Dimension Dim>
class Distributor {
public:
  /**
    Distributor Constructor: Don't use member references!
    Distributor must be constructed before any other members
  **/
  Distributor(const Parameters<Real,Dim>& params,
              Particles<Real,Dim>& particles_):
    environment{false}, comm_active{comm_world, boost::mpi::comm_attach},
    parameters{params}, particles{particles_} {};

  ~Distributor()                             = default;
  Distributor(const Distributor&)            = delete;
  Distributor& operator=(const Distributor&) = delete;
  Distributor(Distributor&&) noexcept        = delete;
  Distributor& operator=(Distributor&&)      = delete;

  /**
    Handle setting initial node boundaries
    required as members not being initilized when constructor is called
  **/
  void Initilize() {
    this->DistributeFluid(parameters.GetInitialFluid());
  }

  /**
    Set node bounds based upon equal spacing in AABB
  **/
  void SetNodeBounds(const AABB<Real,Dim>& aabb) {
    Real node_length = aabb.Length()/comm_active.size();
    node_start = comm_active.rank() * node_length;
    node_end = node_start + node_length;
  }

  /**
    Return rank of node to the left or MPI_PROC_NULL
  **/
  int NodeToLeft() const {
    int left = MPI_PROC_NULL;
    if(comm_active.rank() != 0)
      left = comm_active.rank() - 1;

    return left;
  }

  /**
    Return rank of node to the right or MPI_PROC_NULL
  **/
  int NodeToRight() const {
    int right = MPI_PROC_NULL;
    if(comm_active.rank() != comm_active.size()-1)
      right = comm_active.rank() + 1;

    return right;
  }

  /**
    Construct water volume spread across multiple nodes
  **/
  void DistributeFluid(const AABB<Real, Dim>& global_fluid) {
    this->SetNodeBounds(global_fluid);

    // Create a bounding box with length dimension that coincides with particle spacing
    // Assume particles will be placed at edge of boundary
    AABB<Real,Dim> local_fluid{global_fluid};
    Real node_length = node_end - node_start;
    int particle_count_x = floor(node_length/parameters.GetParticleSpacing());
    if(particle_count_x < 2)
      throw std::string("Not enough particles per node\n");
    Real local_fluid_length = (particle_count_x-1) * parameters.GetParticleSpacing();
    local_fluid.min.x = global_fluid.min.x + local_fluid_length*comm_active.rank();
    local_fluid.max.x = local_fluid.min.x + local_fluid_length;

    // Fill AABB with particles
    particles.AddFluid(local_fluid);
  }

private:
  boost::mpi::environment environment;
  boost::mpi::communicator comm_world;
  const boost::mpi::communicator& comm_active;
  const Parameters<Real,Dim>& parameters;
  Particles<Real,Dim>& particles;
  Real node_start;
  Real node_end;
};
