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
  Real length() const {
    return max.x - min.x;
  }

  /**
    @brief Return y extent of AABB
  **/
  Real height() const {
    return max.y - min.y;
  }

  /**
    @brief return area of 2D AABB
  **/
  Real area() const {
    return this->length() * this->height();
  }

  /**
    @brief return area of 2D AABB
  **/
  Real volume() const {
    return this->area();
  }

  /**
    @brief return length, width vector
  **/
  Vec<Real,two_dimensional> extent() const {
    return max - min;
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
  Real length() const {
    return max.x - min.x;
  }

  /**
    @brief Return y extent of AABB
  **/
  Real height() const {
    return max.y - min.y;
  }

  /**
    @brief Return z extent of AABB
  **/
  Real depth() const {
    return max.z - min.z;
  }

  /**
    @brief return volume of 3D AABB
  **/
  Real volume() const {
    return this->length() * this->height() * this->depth();
  }

  /**
    @brief return length, width, depth vector
  **/
  Vec<Real,three_dimensional> extent() const {
    return max - min;
  }

};

/**
   @brief return a vec containing the number of bins
   that the aabb can be divided into given a given spacing
**/
template <typename Real, Dimension Dim>
Vec<std::size_t, Dim> bin_count_in_volume(const AABB<Real, Dim>& aabb, Real spacing) {
  const Vec<Real,Dim> space_counts{floor((aabb.max - aabb.min) / spacing)};
  Vec<std::size_t,Dim> int_counts(static_cast< Vec<std::size_t,Dim> >(space_counts));
    return int_counts;
  }

