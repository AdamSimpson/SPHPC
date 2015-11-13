#pragma once

#include <boost/mpi.hpp>
#include <utility>

template<typename Real, Dimension Dim>
class Distributor {
public:
  Distributor() : environment{false} {};

  ~Distributor()                             = default;
  Distributor(const Distributor&)            = delete;
  Distributor& operator=(const Distributor&) = delete;
  Distributor(Distributor&&) noexcept        = delete;
  Distributor& operator=(Distributor&&)      = delete;

private:
  boost::mpi::environment environment;
  boost::mpi::communicator comm_world;
};
