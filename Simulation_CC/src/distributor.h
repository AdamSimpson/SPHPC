#pragma once

#include <boost/mpi.hpp>

template<typename Real, typename Integer, Dimension Dim>
class Distributor {
public:
  Distributor()                              = default;
  ~Distributor()                             = default;
  Distributor(const Distributor&)            = default;
  Distributor& operator=(const Distributor&) = default;
  Distributor(Distributor&&) noexcept        = default;
  Distributor& operator=(Distributor&&)      = default;
  
private:
};
