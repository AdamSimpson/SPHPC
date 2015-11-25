#pragma once

template<typename Real, Dimension Dim>
class Particles;

template<typename Real, Dimension Dim>
class Neighbors {
public:
  Neighbors(const Particles<Real, Dim>& particles): m_particles{particles} {};

  ~Neighbors()                           = default;
  Neighbors(const Neighbors&)            = default;
  Neighbors& operator=(const Neighbors&) = default;
  Neighbors(Neighbors&&) noexcept        = default;
  Neighbors& operator=(Neighbors&&)      = default;

private:
  const Particles<Real, Dim>& m_particles;
};
