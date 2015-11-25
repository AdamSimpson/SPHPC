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
    due to MPI_Init() call.

    Requires Initilization step after construction of member references

    compute_comm is split with color value 1 to differentiate from I/O and viz
  **/
  Distributor(const Parameters<Real,Dim>& parameters,
              Particles<Real,Dim>& particles):
    m_environment{false}, m_comm_compute{m_comm_world.split(1)},
    m_parameters{parameters}, m_particles{particles} {};

  ~Distributor()                             = default;
  Distributor(const Distributor&)            = delete;
  Distributor& operator=(const Distributor&) = delete;
  Distributor(Distributor&&) noexcept        = delete;
  Distributor& operator=(Distributor&&)      = delete;

  /**
    @brief Handle initilization requiring constructed member references
  **/
  void Initilize() {
    this->SetNodeBounds(m_parameters.GetInitialFluid());
  }

  /**
    @brief Set node bounds based upon equal spacing in AABB
  **/
  void SetNodeBounds(const AABB<Real,Dim>& aabb) {
    Real node_length = aabb.Length()/m_comm_compute.size();
    m_node_start = m_comm_compute.rank() * node_length;

    if(this->LastNode())
      m_node_end = m_parameters.GetBoundary().max.x;
    else
      m_node_end = m_node_start + node_length;
  }

  /**
    @return Return the rank of the last rank in compute communicator
  **/
  bool LastNode() const {
    return m_comm_compute.rank() == m_comm_compute.size()-1;
  }

  /**
    @return Return rank of node to the left or MPI_PROC_NULL
  **/
  int NodeToLeft() const {
    int left = MPI_PROC_NULL;
    if(m_comm_compute.rank() != 0)
      left = m_comm_compute.rank() - 1;

    return left;
  }

  /**
    @brief Return rank of node to the right or MPI_PROC_NULL
  **/
  int NodeToRight() const {
    int right = MPI_PROC_NULL;
    if(m_comm_compute.rank() != m_comm_compute.size()-1)
      right = m_comm_compute.rank() + 1;

    return right;
  }

  /**
    @brief getter for compute  MPI communicator
  **/
  boost::mpi::communicator GetComputeComm() const {
    return m_comm_compute;
  }

  /**
    @brief getter for compute integer rank
  **/
  int GetComputeRank() const {
    return m_comm_compute.rank();
  }

  /**
    Construct water volume spread across multiple nodes
  **/
  void DistributeFluid(const AABB<Real, Dim>& global_fluid) {
    // Create a bounding box with length dimension that coincides with particle spacing
    // |--|--|--|
    // -o--o--o-
    Real spacing = m_parameters.GetParticleRestSpacing();
    AABB<Real,Dim> local_fluid{global_fluid};
    size_t global_particle_count_x = floor(global_fluid.Length()/spacing);

    // Append particles on last rank if they can't be evenly distributed
    int particle_count_x = global_particle_count_x / m_comm_compute.size();
    if(this->LastNode()) {
      int remaining = global_particle_count_x - (particle_count_x*m_comm_compute.size());
      particle_count_x += remaining;
    }

    if(particle_count_x < 2)
      throw std::runtime_error("Less than two particles in x dimension!");

    Real local_fluid_length = (particle_count_x) * spacing;
    local_fluid.min.x = global_fluid.min.x + local_fluid_length*m_comm_compute.rank();
    // small delta added to ensure no round off error
    local_fluid.max.x = local_fluid.min.x + local_fluid_length + spacing/10.0;

    // Fill AABB with particles
    m_particles.AddFluid(local_fluid);
  }

  std::size_t GetGlobalParticleCount() const {
    std::size_t global_count = 0;
    boost::mpi::all_reduce(m_comm_compute, m_particles.GetLocalCount(), global_count, std::plus<std::size_t>());
    return global_count;
  }

  void ExchangeOOB() {

  }

  void ExchangeHalo() {

  }

private:
  boost::mpi::environment m_environment;
  boost::mpi::communicator m_comm_world;          /**< world communicator **/
  const boost::mpi::communicator m_comm_compute;  /**< compute subset of simulation **/
  const Parameters<Real,Dim>& m_parameters;
  Particles<Real,Dim>& m_particles;
  Real m_node_start;                              /** x coordinate of node domain start **/
  Real m_node_end;                                /** x coordinate of node domain end   **/
};

/**
  Boost Vec<Real,Dim> optimizations
**/

/**
  BOOST_IS_MPI_DATATYPE(Vec<Real,Dim>)
**/
namespace boost { namespace mpi {
  template <typename Real, Dimension Dim>
  struct is_mpi_datatype< Vec<Real,Dim> > : mpl::true_ { };
}}

/**
  BOOST_CLASS_TRACKING(Vec<Real,Dim>,track_never)
**/
namespace boost { namespace serialization {
template <typename Real, Dimension Dim>
struct tracking_level< Vec<Real,Dim> >
{
    typedef mpl::integral_c_tag tag;
    typedef mpl::int_<track_never> type;
    BOOST_STATIC_CONSTANT(
        int,
        value = tracking_level::type::value
    );
}; }}

/**
  BOOST_CLASS_IMPLEMENTATION(Vec<Real,Dim>, object_serializable)
**/
namespace boost { namespace serialization {
template <typename Real, Dimension Dim>
struct implementation_level< Vec<Real,Dim> >
{
    typedef mpl::integral_c_tag tag;
    typedef mpl::int_<object_serializable> type;
    BOOST_STATIC_CONSTANT(
        int,
        value = implementation_level::type::value
    );
}; }}

/**
BOOST_IS_BITWISE_SERIALIZABLE(Vec<Real,Dim>)
**/
namespace boost { namespace serialization {
    template <typename Real, Dimension Dim>
    struct is_bitwise_serializable< Vec<Real,Dim> > : mpl::true_ {};
} }
