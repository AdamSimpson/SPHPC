#pragma once

#include "dimension.h"
#include "vec.h"

/**
  Class to handle Axis aligned boundary box and 2D equivalent
**/

template<typename Real, Dimension dim> class AABB;

template<typename Real>
class AABB<Real, two_dimensional> {
public:
  AABB(): min(0.0), max(0.0) {};

  Vec<Real, two_dimensional> min;
  Vec<Real, two_dimensional> max;

  Real Length() const {
    return max.x - min.x;
  }

  Real Height() const {
    return max.y - min.y;
  }

  Real Volume() const {
    return this->Length() * this->Height();
  }

};

template<typename Real>
class AABB<Real, three_dimensional> {
public:
  AABB(): min(0.0), max(0.0) {};

  Vec<Real, three_dimensional> min;
  Vec<Real, three_dimensional> max;

  Real Length() const {
    return max.x - min.x;
  }

  Real Height() const {
    return max.y - min.y;
  }

  Real Depth() const {
    return max.z - min.z;
  }

  Real Volume() const {
    return this->Length() * this->Height() * this->Depth();
  }
};
