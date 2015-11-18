#pragma once

#include "dimension.h"

template<typename Real, Dimension Dim>
class Kernels {
public:
  Kernels(const Parameters<Real, Dim>& params): parameters(params) {};
  Real Poly6() {};

  private:
    const Parameters<Real,Dim>& parameters;
};
