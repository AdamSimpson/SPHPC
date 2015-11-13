#pragma once

template<typename Real, Dimension Dim>
class Particles;

template<typename Real, Dimension Dim>
class Neighbors {
public:
  Neighbors(const Particles<Real, Dim>& particles_): particles{particles_} {};

  ~Neighbors()                           = default;
  Neighbors(const Neighbors&)            = default;
  Neighbors& operator=(const Neighbors&) = default;
  Neighbors(Neighbors&&) noexcept        = default;
  Neighbors& operator=(Neighbors&&)      = default;

private:
  const Particles<Real, Dim>& particles;
};
