#include <iostream>

/**
A Vector class to handle arbitrary type and dimension
**/

// Determine the alignment of struct - requires C++14
constexpr int vec_alignment(const size_t struct_bytes) {
 int alignments[5] = {1, 2, 3, 8, 16};

 // Get first available alignment >= to struct size
 int first_larger = -1;
 for(int alignment : alignments) {
   first_larger = alignment;
   if(first_larger >= struct_bytes)
     break;
 }

 return first_larger;
}

template <typename T, int n>
struct Vector {
 T data[n];

 constexpr Vector(std::initializer_list<T> l){ for(int i=0; i<l.size(); ++i){ data[i] = l[i];} }
 constexpr explicit Vector(T *data_){ for(int i=0; i<n; ++i){ data[i] = data_[i]; } }
 constexpr explicit Vector(const T value){ for(int i=0; i<n; ++i){ data[i] = value;} }

 T& operator[] (const size_t index) { return data[index]; }
 const T& operator[] (const size_t index) const { return data[index]; }

};

/**
Specialization for 2,3 dimensional arrays
Could use partial specialization on T but would loose alignment without some trick
TO DO : impliment partial specialization
**/
template <typename T>
struct Vector<T, 2> {
 union {
   T data[2];
   struct { T x, y; };
 };

 Vector(const float x_, const float y_) : x(x_), y(y_) {}
 constexpr explicit Vector(const T value) : x(value), y(value) {}
 constexpr explicit Vector(T *data_) : x(data_[0]), y(data_[1]) {}

 T& operator[] (const size_t index) { return data[index]; }
 const T& operator[] (const size_t index) const { return data[index]; }
};

template <typename T>
struct Vector<T, 3> {
 union {
   T data[3];
   struct { T x, y, z; };
   struct { T r, g, b; };
 };

 Vector(const float x_, const float y_, const float z_) : x(x_), y(y_), z(z_) {}
 constexpr explicit Vector(const T value) : x(value), y(value), z(value) {}
 constexpr explicit Vector(T *data_) : x(data_[0]), y(data_[1]), z(data_[2]) {}

 T& operator[] (const size_t index) { return data[index]; }
 const T& operator[] (const size_t index) const { return data[index]; }
};

/**
Free Operators
**/

template<typename T, int n>
Vector<T,n> operator+(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 T data[n];
 for(int i=0; i<n; i++)
   data[i] = lhs[i] + rhs[i];

 return Vector<T,n>(data);
}

template<typename T, int n>
Vector<T,n> operator-(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 T data[n];
 for(int i=0; i<n; i++)
   data[i] = lhs[i] - rhs[i];

 return Vector<T,n>(data);
}

template<typename T, int n>
Vector<T,n> operator*(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 T data[n];
 for(int i=0; i<n; i++)
   data[i] = lhs[i] * rhs[i];

 return Vector<T,n>(data);
}

template<typename T, int n>
Vector<T,n> operator*(const Vector<T,n>& lhs, const T rhs) {
 T data[n];
 for(int i=0; i<n; i++)
   data[i] = lhs[i] * rhs;

 return Vector<T,n>(data);
}

template<typename T, int n>
Vector<T,n> operator*(const T lhs, const Vector<T,n>& rhs) {
 T data[n];
 for(int i=0; i<n; i++)
   data[i] = lhs * rhs[i];

 return Vector<T,n>(data);
}

template<typename T, int n>
Vector<T,n>& operator+=(Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] += rhs[i];

 return lhs;
}

template<typename T, int n>
Vector<T,n>& operator-=(Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] -= rhs[i];

 return lhs;
}

template<typename T, int n>
Vector<T,n>& operator*=(Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] *= rhs[i];

 return lhs;
}

template<typename T, int n>
Vector<T,n>& operator*=(Vector<T,n>& lhs, const T& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] *= rhs;

 return lhs;
}

template<typename T, int n>
Vector<T,n>& operator/=(Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] /= rhs[i];

 return lhs;
}

template<typename T, int n>
Vector<T,n>& operator/=(Vector<T,n>& lhs, const T& rhs) {
 for(int i=0; i<n; i++)
   lhs[i] /= rhs;

 return lhs;
}

template<typename T, int n>
T Dot(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 T dot;
 for(int i=0; i<n; ++i)
   dot += lhs[i] * rhs[i];

   return dot;
}

template<typename T, int n>
T MagnitudeSquared(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 return Dot(lhs, rhs);
}

template<typename T, int n>
T Magnitude(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 return sqrt(MagnitudeSquared(lhs, rhs));
}

template<typename T, int n>
T InverseMagnitude(const Vector<T,n>& lhs, const Vector<T,n>& rhs) {
 return static_cast<T>(1.0) / sqrt(MagSquared(lhs, rhs));
}
