#pragma once

#include "dimension.h"
#include <string>

template<typename Real, typename Integer, Dimension Dim>
class Parameters {
public:
  Parameters(const std::string& file_name) {
    max_particles_local = 1000;
    printf("params and shit\n");
  };

  ~Parameters()                            = default;
  Parameters(const Parameters&)            = default;
  Parameters& operator=(const Parameters&) = default;
  Parameters(Parameters&&) noexcept        = default;
  Parameters& operator=(Parameters&&)      = default;

  Integer GetMaxParticlesLocal() const {
    return max_particles_local;
  };

private:
  Integer max_particles_local;
};
