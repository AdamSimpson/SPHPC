#pragma once

#include <iostream>
#include "utility_math.h"
#include "dimension.h"
#include "hemi/hemi.h"

// Forward decleration of Axis Aligned Boundary Box
template<typename Real, Dimension dim> class AABB;

/////////////////////////////////////////////////
// A basic vector class to handle arbitrary type
// and dimensions
/////////////////////////////////////////////////

/**
  @brief Determine the alignment of struct - requires C++14 to work correctly
**/
constexpr int vec_alignment(const size_t struct_bytes) {
/*
  // Available CUDA alignments
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
  return 8;
}

/**
  @struct Vec
  @brief Generic n-dim Vec implementation
  @tparam T component type
  @tparam n dimension
**/
template <typename T, int N>
struct Vec {
  /**
    @brief Component storage
  **/
  T m_data[N];

  /**
    @brief Default constructor: components uninitialized
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  Vec() = default;

  /**
    @brief Initilizer list constructor
    @param[in] l initilizer list
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  explicit Vec(std::initializer_list<T> l){ for(int i=0; i<l.size() && i<N; ++i){ m_data[i] = l[i];} }

  /**
    @brief Pointer constructor: components initilized to pointed data
    @param[in] p_data pointer to data
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  explicit Vec(const T *const p_data){ for(int i=0; i<N; ++i){ m_data[i] = p_data[i]; } }

  /**
    @brief Scalar constructor: components initilized to scalar value
    @param[in] value scalar component value
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  explicit Vec(const T value){ for(int i=0; i<N; ++i){ m_data[i] = value;} }

  /**
    @brief Cast operator: static_cast() components of n-dim Vec
  **/
  template<typename T_out>
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  operator Vec<T_out,N>() {
    T_out data_out[N];
    for(int i=0; i<N; ++i) {
      data_out[i] = static_cast<T_out>(m_data[i]);
    }
    return Vec<T_out,N>(data_out);
  }

  /**
    @brief Subscript operator: access vector components with bracket notation
    @param[in] index of element
    @return index'th element of Vec
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  T& operator[] (const size_t index) { return m_data[index]; }

  /**
    @brief Const subscript operator: access vector components with bracket notation
    @param[in] index of element
    @return index'th element of Vec
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  const T& operator[] (const size_t index) const { return m_data[index]; }

};

/**
  @struct Vec
  @brief Specialization for 2-dim Vec
**/
template <typename T>
struct Vec<T, 2> {
  /**
    @brief union component storage
  **/
  union {
    T m_data[2];          /**< component array */
    struct { T x, y; }; /**<  x, y access */
  };

  /**
    @brief Default constructor: components uninitialized
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  Vec() = default;

  /**
    @brief Constructor: components initilized
    @param[in] x_ x component
    @param[in] y_ y component
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  explicit Vec(const T x_, const T y_) : x(x_), y(y_) {}

  /**
    @brief Constructor: components initilized to value
    @param[in] value x, y component
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  constexpr explicit Vec(const T value) : x(value), y(value) {}

  /**
    @brief Pointer constructor: components initilized to pointed data
    @param[in] p_data pointer to data
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  constexpr explicit Vec(T *data_in) : x(data_in[0]), y(data_in[1]) {}

  /**
    @brief Construct Vec2 from Vec3 by truncating z
    @param[in] vec input Vec3
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  constexpr explicit Vec(const Vec<T,3>& vec): x{vec.x}, y{vec.y} {}

  /**
    @brief Cast operator: static_cast() data of Vec2
  **/
  template<typename T_out>
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  operator Vec<T_out,2>() const {
    return Vec<T_out,2>(static_cast<T_out>(x), static_cast<T_out>(y)); }

  /**
    @brief Subscript operator: access vector components with bracket notation
    @param[in] index of element
    @return index'th element of Vec
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  T& operator[] (const size_t index) { return m_data[index]; }

  /**
    @brief Const subscript operator: access vector components with bracket notation
    @param[in] index of element
    @return index'th element of Vec
  **/
  HEMI_DEV_CALLABLE_INLINE_MEMBER
  const T& operator[] (const size_t index) const { return m_data[index]; }
};

/**
  @struct Vec
  @brief Specialization for 3-dim Vec
**/
template <typename T>
struct Vec<T, 3> {

  union {
    T m_data[3];              /**< component array */
    struct { T x, y, z; };  /**< x, y, z access  */
    struct { T r, g, b; };  /**< r, g, b access  */
  };

  /**
    @brief Default constructor: components uninitialized
  **/
  Vec() = default;

  /**
    @brief Constructor: components initilized
    @param[in] x_ x component
    @param[in] y_ y component
    @param[in] z_ z component
  **/
  explicit Vec(const T x_, const T y_, const T z_) : x(x_), y(y_), z(z_) {}

  /**
    @brief Constructor: components initilized to value
    @param[in] value x, y, z  component
  **/
  constexpr explicit Vec(const T value) : x(value), y(value), z(value) {}

  /**
    @brief Pointer constructor: components initilized to pointed data
    @param[in] p_data pointer to data
  **/
  constexpr explicit Vec(T *data_) : x(data_[0]), y(data_[1]), z(data_[2]) {}

  /**
    @brief Construct Vec2 from Vec3 by setting z to 0
    @param[in] vec Vec2
  **/
  explicit Vec(const Vec<T, 2>& vec_): x{vec_[0]}, y{vec_[1]}, z{0} {}

  /**
    @brief Construct Vec2 from Vec3 by specifying z
    @param[in] vec Vec2 x,y values
    @param[in] z z value
  **/
  explicit Vec(const Vec<T, 2>& vec_, T z_): x{vec_[0]}, y{vec_[1]}, z{z_} {}

  /**
    @brief Cast operator: static_cast() data of Vec3
  **/
  template<typename T_out> operator Vec<T_out,3>() const
    { return Vec<T_out,3>(static_cast<T_out>(x), static_cast<T_out>(y), static_cast<T_out>(z)); }

  /**
    @brief Subscript operator: access vector components with bracket notation
    @param[in] index of element
    @return index'th element of Vec
  **/
  T& operator[] (const size_t index) { return m_data[index]; }

  /**
    @brief Const subscript operator: access vector components with bracket notation
    @param[in] index of element
    @return index'th element of Vec
  **/
  const T& operator[] (const size_t index) const { return m_data[index]; }
};

/**
Free Operators
**/

/**
  @brief vec vec component wise addition operator
**/
template<typename T, int n>
Vec<T,n> operator+(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] + rhs[i];

  return Vec<T,n>(data);
}

/**
  @brief vec scalar component wise addition operator
**/
template<typename T, int n>
Vec<T,n> operator+(const Vec<T,n>& lhs, const T rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] + rhs;

  return Vec<T,n>(data);
}

/**
  @brief vec vec component wise subtraction operator
**/
template<typename T, int n>
Vec<T,n> operator-(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] - rhs[i];

  return Vec<T,n>(data);
}

/**
  @brief vec scalar component wise subtraction operator
**/
template<typename T, int n>
Vec<T,n> operator-(const Vec<T,n>& lhs, const T rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] - rhs;

  return Vec<T,n>(data);
}

/**
  @brief vec vec component wise multiplication operator
**/
template<typename T, int n>
Vec<T,n> operator*(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] * rhs[i];

  return Vec<T,n>(data);
}

/**
  @brief vec scalar component wise multiplication operator
**/
template<typename T, int n>
Vec<T,n> operator*(const Vec<T,n>& lhs, const T rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] * rhs;

  return Vec<T,n>(data);
}

/**
  @brief scalar vec component wise multiplication operator
**/
template<typename T, int n>
Vec<T,n> operator*(const T lhs, const Vec<T,n>& rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs * rhs[i];

  return Vec<T,n>(data);
}

/**
  @brief vec scalar component wise divison operator
**/
template<typename T, int n>
Vec<T,n> operator/(const Vec<T,n>& lhs, const T rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] / rhs;

  return Vec<T,n>(data);
}

/**
  @brief vec vec component wise accumulate operator
**/
template<typename T, int n>
Vec<T,n>& operator+=(Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] += rhs[i];

  return lhs;
}

/**
  @brief vec scalar component wise accumulate operator
**/
template<typename T, int n>
Vec<T,n>& operator+=(Vec<T,n>& lhs, const T rhs) {
  for(int i=0; i<n; i++)
    lhs[i] += rhs;

  return lhs;
}

/**
  @brief vec vec component wise decrement operator
**/
template<typename T, int n>
Vec<T,n>& operator-=(Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] -= rhs[i];

  return lhs;
}

/**
  @brief vec scalar component wise decrement operator
**/
template<typename T, int n>
Vec<T,n>& operator-=(Vec<T,n>& lhs, const T rhs) {
  for(int i=0; i<n; i++)
    lhs[i] -= rhs;

  return lhs;
}

/**
  @brief vec vec component wise multiplication assignment operator
**/
template<typename T, int n>
Vec<T,n>& operator*=(Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] *= rhs[i];

  return lhs;
}

/**
  @brief vec vec component wise multiplication assignment operator
**/
template<typename T, int n>
Vec<T,n>& operator*=(Vec<T,n>& lhs, const T& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] *= rhs;

  return lhs;
}

/**
  @brief vec vec component wise division assignment operator
**/
template<typename T, int n>
Vec<T,n>& operator/=(Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] /= rhs[i];

  return lhs;
}

/**
  @brief vec scalar component wise division assignment operator
**/
template<typename T, int n>
Vec<T,n>& operator/=(Vec<T,n>& lhs, const T& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] /= rhs;

  return lhs;
}

/**
  @brief vector dot product
**/
template<typename T, int n>
T Dot(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  T dot;
  for(int i=0; i<n; ++i)
    dot += lhs[i] * rhs[i];

  return dot;
}

/**
  @brief vector magnitude squared
**/
template<typename T, int n>
T MagnitudeSquared(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  return Dot(lhs, rhs);
}

/**
  @brief vector magnitude
**/
template<typename T, int n>
T Magnitude(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  return sqrt(MagnitudeSquared(lhs, rhs));
}

/**
  @brief reciprical vector magnitude
**/
template<typename T, int n>
T InverseMagnitude(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  return static_cast<T>(1.0) / sqrt(MagSquared(lhs, rhs));
}

/**
  @brief component wise floor
**/
template<typename T, int n>
Vec<T,n> Floor(const Vec<T,n>& vec) {
   Vec<T,n> floor_vec;
   for(int i=0; i<n; i++)
     floor_vec[i] = floor(vec[i]);

  return floor_vec;
}

/**
  @brief component wise vector sum
**/
template<typename T, int n>
T Sum(const Vec<T,n>& vec) {
  T sum = 0;
  for(int i=0; i<n; ++i)
    sum += vec[i];

  return sum;
}

/**
  @brief component wise clamp
  @return newly constructed clamped vector
**/
template<typename T, int n>
Vec<T,n> Clamp(const Vec<T,n>& vec, const Vec<T,n>& lower, const Vec<T,n>& upper) {
   Vec<T,n> clamp_vec;
   for(int i=0; i<n; i++)
     clamp_vec[i] = Utility::Clamp(vec[i], lower[i], upper[i]);

  return clamp_vec;
}

/**
  @brief component wise in-place clamp
  @return original vector with values clamped
**/
template<typename T, int n>
void ClampInPlace(Vec<T,n>& vec, const Vec<T,n>& lower, const Vec<T,n>& upper) {
   for(int i=0; i<n; i++)
     Utility::ClampInPlace(vec[i], lower[i], upper[i]);
}
