#pragma once

template<typename Real, typename Integer, Dimension Dim>
class Particles;

template<typename Real, typename Integer, Dimension Dim>
class Neighbors {
public:
  Neighbors(const Particles<Real, Integer, Dim>& particles_): particles{particles_} {};

  ~Neighbors()                           = default;
  Neighbors(const Neighbors&)            = default;
  Neighbors& operator=(const Neighbors&) = default;
  Neighbors(Neighbors&&) noexcept        = default;
  Neighbors& operator=(Neighbors&&)      = default;

private:
  const Particles<Real, Integer, Dim>& particles;
};
