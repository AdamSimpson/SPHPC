#ifndef SPH_SRC_VEC_H_
#define SPH_SRC_VEC_H_

#include <iostream>
#include "hemi/hemi.h"

/**
A Vector class to handle arbitrary type and dimension
**/

// Determine the alignment of struct - requires C++14

/* C++14 */
constexpr int vec_alignment(/*const size_t struct_bytes*/) {
/*
 int alignments[5] = {1, 2, 3, 8, 16};

 // Get first available alignment >= to struct size
 int first_larger = -1;
 for(int alignment : alignments) {
   first_larger = alignment;
   if(first_larger >= struct_bytes)
     break;
 }

 return first_larger;
 */
  return 16;
}

template <typename T, int n>
struct HEMI_ALIGN(vec_alignment(/*sizeof(T)*n*/)) Vector {
 T data[n];

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 constexpr Vector() {};

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 /*constexpr C++14*/ Vector(std::initializer_list<T> l){ for(int i=0; i<l.size(); ++i){ data[i] = l[i];} }

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 /*constexpr C++14*/explicit Vector(T *data_){ for(int i=0; i<n; ++i){ data[i] = data_[i]; } }

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 /*constexpr C++14*/explicit Vector(const T value){ for(int i=0; i<n; ++i){ data[i] = value;} }

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 T& operator[] (const size_t index) { return data[index]; }

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 const T& operator[] (const size_t index) const { return data[index]; }
};

/**
Specialization for 2,3 dimensional arrays
Could use partial specialization on T but would loose alignment without some trick
TO DO : impliment partial specialization
**/
template <typename T>
struct HEMI_ALIGN(vec_alignment(/*sizeof(T)*2*/)) Vector<T, 2> {
 union {
   T data[2];
   struct { T x, y; };
 };

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 /*constexpr NVCC-broke*/ Vector() {};

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 constexpr Vector(const float x_, const float y_) : x(x_), y(y_) {}

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 constexpr explicit Vector(const T value) : x(value), y(value) {}

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 constexpr explicit Vector(T *data_) : x(data_[0]), y(data_[1]) {}

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 T& operator[] (const size_t index) { return data[index]; }

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 const T& operator[] (const size_t index) const { return data[index]; }
};

template <typename T>
struct HEMI_ALIGN(vec_alignment(/*sizeof(T)*3*/)) Vector<T, 3> {
 union {
   T data[3];
   struct { T x, y, z; };
   struct { T r, g, b; };
 };

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 /*constexpr NVCC-broke*/ Vector() {};

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 constexpr Vector(const float x_, const float y_, const float z_) : x(x_), y(y_), z(z_) {}

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 constexpr explicit Vector(const T value) : x(value), y(value), z(value) {}

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 constexpr explicit Vector(T *data_) : x(data_[0]), y(data_[1]), z(data_[2]) {}

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 T& operator[] (const size_t index) { return data[index]; }

 HEMI_DEV_CALLABLE_INLINE_MEMBER
 const T& operator[] (const size_t index) const { return data[index]; }
};

/**
Free Operators
**/

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
Vector<T,n> operator+(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 T data[n];
 for(int i=0; i<n; i++)
   data[i] = lhs[i] + rhs[i];

 return Vector<T,n>(data);
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
Vector<T,n> operator-(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 T data[n];
 for(int i=0; i<n; i++)
   data[i] = lhs[i] - rhs[i];

 return Vector<T,n>(data);
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
Vector<T,n> operator*(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 T data[n];
 for(int i=0; i<n; i++)
   data[i] = lhs[i] * rhs[i];

 return Vector<T,n>(data);
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
Vector<T,n> operator*(const Vector<T,n>& lhs, const T rhs) {
 T data[n];
 for(int i=0; i<n; i++)
   data[i] = lhs[i] * rhs;

 return Vector<T,n>(data);
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
Vector<T,n> operator*(const T lhs, const Vector<T,n>& rhs) {
 T data[n];
 for(int i=0; i<n; i++)
   data[i] = lhs * rhs[i];

 return Vector<T,n>(data);
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
Vector<T,n>& operator+=(Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] += rhs[i];

 return lhs;
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
Vector<T,n>& operator-=(Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] -= rhs[i];

 return lhs;
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
Vector<T,n>& operator*=(Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] *= rhs[i];

 return lhs;
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
Vector<T,n>& operator*=(Vector<T,n>& lhs, const T& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] *= rhs;

 return lhs;
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
Vector<T,n>& operator/=(Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] /= rhs[i];

 return lhs;
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
Vector<T,n>& operator/=(Vector<T,n>& lhs, const T& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] /= rhs;

 return lhs;
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
T Dot(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 T dot;
 for(int i=0; i<n; ++i)
   dot += lhs[i] * rhs[i];

   return dot;
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
T MagnitudeSquared(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 return Dot(lhs, rhs);
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
T Magnitude(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 return sqrt(MagnitudeSquared(lhs, rhs));
}

template<typename T, int n>
HEMI_DEV_CALLABLE_INLINE
T InverseMagnitude(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 return static_cast<T>(1.0) / sqrt(MagSquared(lhs, rhs));
}

#endif
