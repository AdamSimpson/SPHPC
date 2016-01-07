#pragma once

#include "dimension.h"
#include "vec.h"
#include <math.h>

/**
   Empty generic template class
 **/
template<typename Real, Dimension Dim> class Poly6;
template<typename Real, Dimension Dim> class Del_Poly6;
template<typename Real, Dimension Dim> class Del_Spikey;
template<typename Real, Dimension Dim> class C_Spline;

/**
   3D poly6 specilization
 **/
template <typename Real>
class Poly6<Real, three_dimensional> {
public:
  Poly6(const Real smoothing_radius): h_{smoothing_radius},
                                      norm_{static_cast<Real>(315.0/(64.0*M_PI*pow(h_, 9.0)))} {};

  Real operator()(const Real r) const {
    if(r > h_)
      return 0.0;

    return norm_ * (h_*h_ - r*r) * (h_*h_ - r*r) * (h_*h_ - r*r);
  };

private:
  Real h_;
  Real norm_;
};

/**
   2D poly6 specilization
**/
template <typename Real>
class Poly6<Real, two_dimensional> {
public:
  Poly6(const Real smoothing_radius): h_{smoothing_radius},
                                      norm_{static_cast<Real>(4.0/(M_PI*pow(h_, 8.0)))} {};

  Real operator()(const Real r_mag) const {
    if(r_mag > h_)
     return 0.0;

    return norm_ * (h_*h_ - r_mag*r_mag) * (h_*h_ - r_mag*r_mag) * (h_*h_ - r_mag*r_mag);
  };

private:
  Real h_;
  Real norm_;
};

/**
   3D gradient poly6 specilization
**/
template <typename Real>
class Del_Poly6<Real, three_dimensional> {
public:
  Del_Poly6(const Real smoothing_radius): h_{smoothing_radius},
                                          norm_{static_cast<Real>(-3.0/2.0 * 315.0/(64.0*M_PI*pow(h_, 9.0)))} {};

  Vec<Real,three_dimensional> operator()(const Vec<Real,three_dimensional>& vec_p,
                                         const Vec<Real,three_dimensional>& vec_q) const {
    const Vec<Real,three_dimensional> r = vec_p - vec_q;
    const Real r_mag_squared = magnitude_squared(r);
    const Real h_squared = h_*h_;
    if(r_mag_squared > h_squared)
      return 0.0;

    // r_mag / r_mag cancels out and is not included
    return norm_ * (h_squared - r_mag_squared) * (h_squared - r_mag_squared) *  r;
  };

private:
  Real h_;
  Real norm_;
};

/**
   2D gradient poly6 specilization
**/
template <typename Real>
class Del_Poly6<Real, two_dimensional> {
public:
  Del_Poly6(const Real smoothing_radius): h_{smoothing_radius},
                                          norm_{static_cast<Real>(-24.0/(M_PI*pow(h_, 8.0)))} {};

  Vec<Real,two_dimensional> operator()(const Vec<Real,two_dimensional>& vec_p,
                                         const Vec<Real,two_dimensional>& vec_q) const {
    const Vec<Real,two_dimensional> r = vec_p - vec_q;
    const Real r_mag_squared = magnitude_squared(r);
    const Real h_squared = h_*h_;

    if(r_mag_squared > h_squared)
      return Vec<Real,two_dimensional>(0.0);

    // r_mag / r_mag cancels out and is not included
    return norm_ * (h_squared - r_mag_squared) * (h_squared - r_mag_squared) *  r;
  };

private:
  Real h_;
  Real norm_;
};

/**
   3D gradient spikey specilization
**/
template <typename Real>
class Del_Spikey<Real, three_dimensional> {
public:
  Del_Spikey(const Real smoothing_radius): h_{smoothing_radius},
                                           norm_{static_cast<Real>(-45.0/(M_PI*pow(h_, 6.0)))} {};

  Vec<Real,three_dimensional> operator()(const Vec<Real,three_dimensional>& vec_p,
                                         const Vec<Real,three_dimensional>& vec_q) const {
    const Vec<Real,three_dimensional> r = vec_p - vec_q;
    const Real r_mag = magnitude(r);
    if(r_mag > h_)
      return Vec<Real,three_dimensional>{0.0};

    return norm_ * (h_ - r_mag) * (h_ - r_mag) *  r / r_mag;
  };

private:
  Real h_;
  Real norm_;
};

/**
   2D gradient spikey specilization
**/
template <typename Real>
class Del_Spikey<Real, two_dimensional> {
public:
  Del_Spikey(const Real smoothing_radius): h_{smoothing_radius},
                                           norm_{static_cast<Real>(-30.0/(M_PI*pow(h_, 5.0)))} {};

  Vec<Real,two_dimensional> operator()(const Vec<Real,two_dimensional>& vec_p,
                                       const Vec<Real,two_dimensional>& vec_q) const {
    const Vec<Real,two_dimensional> r = vec_p - vec_q;
    Real r_mag = magnitude(r);
    if(r_mag > h_)
      return Vec<Real, two_dimensional>{0.0};

    if(r_mag < 0.0001)
      r_mag = 0.0001;

    return norm_ * (h_ - r_mag) * (h_ - r_mag) *  r / r_mag;
  };

private:
  Real h_;
  Real norm_;
};

/**
   3D C_Spline specilization
   norm_ is artificial
**/
template <typename Real>
class C_Spline<Real, three_dimensional> {
public:
  C_Spline(const Real smoothing_radius): h_{smoothing_radius},
                                      norm_{static_cast<Real>(32.0/(M_PI*pow(h_, 9.0)))} {};

  Real operator()(const Real r) const {
    if(r > h_)
      return 0.0;
    else if(r <= h_*static_cast<Real>(0.5))
      return norm_ * (2.0*(h_-r)*(h_-r)*(h_-r)*r*r*r) - (h_*h_*h_*h_*h_*h_/static_cast<Real>(64.0));
    else
      return norm_ * (h_-r)*(h_-r)*(h_-r)*r*r*r;
  };

private:
  Real h_;
  Real norm_;
};

/**
   2D C_Spline specilization
   norm_ is artificial
**/
template <typename Real>
class C_Spline<Real, two_dimensional> {
public:
  C_Spline(const Real smoothing_radius): h_{smoothing_radius},
                                         norm_{static_cast<Real>(32.0/(M_PI*pow(h_, 9.0)))} {};

  Real operator()(const Real r) const {
    if(r > h_)
      return 0.0;
    else if(r <= h_*static_cast<Real>(0.5))
      return norm_ * (2.0*(h_-r)*(h_-r)*(h_-r)*r*r*r) - (h_*h_*h_*h_*h_*h_/static_cast<Real>(64.0));
    else
      return norm_ * (h_-r)*(h_-r)*(h_-r)*r*r*r;
  };

private:
  Real h_;
  Real norm_;
};
