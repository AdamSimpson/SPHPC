#pragma once

#include <iostream>

/**
A Vec class to handle arbitrary type and dimension
**/

// Determine the alignment of struct - requires C++14
constexpr int vec_alignment(const size_t struct_bytes) {
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
  return 8;
}

template <typename T, int n>
struct Vec {
  T data[n];

  Vec() = default;
  /*constexpr C++14*/ explicit Vec(std::initializer_list<T> l){ for(int i=0; i<l.size(); ++i){ data[i] = l[i];} }
  /*constexpr C++14*/ explicit Vec(T *data_){ for(int i=0; i<n; ++i){ data[i] = data_[i]; } }
  /*constexpr C++14*/ explicit Vec(const T value){ for(int i=0; i<n; ++i){ data[i] = value;} }

  /**
    Construct Vec from another Vec, pad with 0's
  **/
/*
  template<typename T_in, int dim_in> Vec(const Vec<T_in, dim_in>& vec_in){
    for(int i=0; i<n; ++i) {
      if(i < dim_in)
        data[i] = static_cast<T>(vec_in[i]);
      else
        data[i] = static_cast<T>(0);
    }
  }
*/

  template<typename T_out> operator Vec<T_out,n>(){
    T_out data_out[n];
    for(int i=0; i<n; ++i) {
      data_out[n] = static_cast<T_out>(data[i]);
    }
    return Vec<T_out,n>(data_out);
  }

  T& operator[] (const size_t index) { return data[index]; }
  const T& operator[] (const size_t index) const { return data[index]; }

};

/**
Specialization for 2,3 dimensional arrays
Could use partial specialization on T but would loose alignment without some trick
TO DO : impliment partial specialization
**/
template <typename T>
struct Vec<T, 2> {
  union {
    T data[2];
    struct { T x, y; };
  };

  Vec() = default;
  explicit Vec(const T x_, const T y_) : x(x_), y(y_) {}
  constexpr explicit Vec(const T value) : x(value), y(value) {}
  constexpr explicit Vec(T *data_) : x(data_[0]), y(data_[1]) {}

  // Copy from Vec3 to Vec2 by truncating z
  constexpr explicit Vec(const Vec<T,3>& vec): x{vec.x}, y{vec.y} {}

  /**
    Construct Vec from another Vec, pad with 0's
  **/
/*
  template<typename T_in, int dim_in> Vec(const Vec<T_in, dim_in>& vec_in){
    for(int i=0; i<2; ++i) {
      if(i < dim_in)
        data[i] = static_cast<T>(vec_in[i]);
      else
        data[i] = static_cast<T>(0);
    }
  }
*/

  template<typename T_out> operator Vec<T_out,2>() const {
    return Vec<T_out,2>(static_cast<T_out>(x), static_cast<T_out>(y)); }

  T& operator[] (const size_t index) { return data[index]; }
  const T& operator[] (const size_t index) const { return data[index]; }
};

template <typename T>
struct Vec<T, 3> {
  union {
    T data[3];
    struct { T x, y, z; };
    struct { T r, g, b; };
  };

  Vec() = default;
  explicit Vec(const T x_, const T y_, const T z_) : x(x_), y(y_), z(z_) {}
  constexpr explicit Vec(const T value) : x(value), y(value), z(value) {}
  constexpr explicit Vec(T *data_) : x(data_[0]), y(data_[1]), z(data_[2]) {}

  /**
    Construct Vec<T,3> from Vec<T,2>
  **/
  explicit Vec(const Vec<T, 2>& vec_): x{vec_[0]}, y{vec_[1]}, z{0} {}
  explicit Vec(const Vec<T, 2>& vec_, T z_): x{vec_[0]}, y{vec_[1]}, z{z_} {}

  template<typename T_out> operator Vec<T_out,3>() const
    { return Vec<T_out,3>(static_cast<T_out>(x), static_cast<T_out>(y), static_cast<T_out>(z)); }

  T& operator[] (const size_t index) { return data[index]; }
  const T& operator[] (const size_t index) const { return data[index]; }
};

/**
Free Operators
**/

template<typename T, int n>
Vec<T,n> operator+(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] + rhs[i];

  return Vec<T,n>(data);
}

template<typename T, int n>
Vec<T,n> operator+(const Vec<T,n>& lhs, const T rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] + rhs;

  return Vec<T,n>(data);
}

template<typename T, int n>
Vec<T,n> operator-(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] - rhs[i];

  return Vec<T,n>(data);
}

template<typename T, int n>
Vec<T,n> operator-(const Vec<T,n>& lhs, const T rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] - rhs;

  return Vec<T,n>(data);
}

template<typename T, int n>
Vec<T,n> operator*(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] * rhs[i];

  return Vec<T,n>(data);
}

template<typename T, int n>
Vec<T,n> operator*(const Vec<T,n>& lhs, const T rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] * rhs;

  return Vec<T,n>(data);
}

template<typename T, int n>
Vec<T,n> operator*(const T lhs, const Vec<T,n>& rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs * rhs[i];

  return Vec<T,n>(data);
}

template<typename T, int n>
Vec<T,n> operator/(const Vec<T,n>& lhs, const T rhs) {
  T data[n];
  for(int i=0; i<n; i++)
    data[i] = lhs[i] / rhs;

  return Vec<T,n>(data);
}

template<typename T, int n>
Vec<T,n>& operator+=(Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] += rhs[i];

  return lhs;
}

template<typename T, int n>
Vec<T,n>& operator+=(Vec<T,n>& lhs, const T rhs) {
  for(int i=0; i<n; i++)
    lhs[i] += rhs;

  return lhs;
}

template<typename T, int n>
Vec<T,n>& operator-=(Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] -= rhs[i];

  return lhs;
}

template<typename T, int n>
Vec<T,n>& operator-=(Vec<T,n>& lhs, const T rhs) {
  for(int i=0; i<n; i++)
    lhs[i] -= rhs;

  return lhs;
}

template<typename T, int n>
Vec<T,n>& operator*=(Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] *= rhs[i];

  return lhs;
}

template<typename T, int n>
Vec<T,n>& operator*=(Vec<T,n>& lhs, const T& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] *= rhs;

  return lhs;
}

template<typename T, int n>
Vec<T,n>& operator/=(Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] /= rhs[i];

  return lhs;
}

template<typename T, int n>
Vec<T,n>& operator/=(Vec<T,n>& lhs, const T& rhs) {
  for(int i=0; i<n; i++)
    lhs[i] /= rhs;

  return lhs;
}

template<typename T, int n>
T Dot(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  T dot;
  for(int i=0; i<n; ++i)
    dot += lhs[i] * rhs[i];

  return dot;
}

template<typename T, int n>
T MagnitudeSquared(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  return Dot(lhs, rhs);
}

template<typename T, int n>
T Magnitude(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  return sqrt(MagnitudeSquared(lhs, rhs));
}

template<typename T, int n>
T InverseMagnitude(const Vec<T,n>& lhs, const Vec<T,n>& rhs) {
  return static_cast<T>(1.0) / sqrt(MagSquared(lhs, rhs));
}

template<typename T, int n>
Vec<T,n> Floor(const Vec<T,n>& vec) {
   Vec<T,n> floor_vec;
   for(int i=0; i<n; i++)
     floor_vec[i] = floor(vec[i]);

  return floor_vec;
}

template<typename T, int n>
T Sum(const Vec<T,n>& vec) {
  T sum = 0;
  for(int i=0; i<n; ++i)
    sum += vec[i];

  return sum;
}
