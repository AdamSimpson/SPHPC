#pragma once

#include "dimension.h"
#include "vec.h"


/////////////////////////////////////////////////
// A basic 2D and 3D Axis aligned boundary box
/////////////////////////////////////////////////

/**
  Empty generic template class
**/
template<typename Real, Dimension dim> class AABB;

/**
  2D AABB specilization
**/
template<typename Real>
class AABB<Real, two_dimensional> {
public:
  /**
    @brief Constructor: initilized components to {0}
  **/
  AABB(): min(0.0), max(0.0) {};

  Vec<Real, two_dimensional> min; /**< minimum x,y,{z} coordinate **/
  Vec<Real, two_dimensional> max; /**< maximum x,y,{z} coordinate **/

  /**
    @brief Return x extent of AABB
  **/
  Real Length() const {
    return max.x - min.x;
  }

  /**
    @brief Return y extent of AABB
  **/
  Real Height() const {
    return max.y - min.y;
  }

  /**
    @brief return area of 2D AABB
  **/
  Real Area() const {
    return this->Length() * this->Height();
  }

  /**
    @brief return area of 2D AABB
  **/
  Real Volume() const {
    return this->Area();
  }

};

/**
  3D AABB specilization
**/
template<typename Real>
class AABB<Real, three_dimensional> {
public:
  /**
    @brief Constructor: initilized components to {0}
  **/
  AABB(): min(0.0), max(0.0) {};

  Vec<Real, three_dimensional> min; /**< minimum x,y,{z} coordinate **/
  Vec<Real, three_dimensional> max; /**< maximum x,y,{z} coordinate **/

  /**
    @brief Return x extent of AABB
  **/
  Real Length() const {
    return max.x - min.x;
  }

  /**
    @brief Return y extent of AABB
  **/
  Real Height() const {
    return max.y - min.y;
  }

  /**
    @brief Return z extent of AABB
  **/
  Real Depth() const {
    return max.z - min.z;
  }

  /**
    @brief return area of 2D AABB
  **/
  Real Volume() const {
    return this->Length() * this->Height() * this->Depth();
  }
};
