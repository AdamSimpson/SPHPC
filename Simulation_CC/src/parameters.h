#pragma once

#include "dimension.h"
#include <string>

template<typename Real, typename Integer, Dimension Dim>
class Parameters {
public:
  Parameters(const std::string& file_name);
  Integer GetMaxParticlesLocal() const { return max_particles_local; };
private:
  Integer max_particles_local;
};

#include "parameters.inl"
