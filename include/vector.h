#ifndef Vec3_H
#define Vec3_H

// This is a placeholder until we replace it with something perhaps from Eigen


template <typename T>
struct Vec3 {
  T x, y, z;

  // Default constructor
  KOKKOS_INLINE_FUNCTION 
  Vec3() : x(0), y(0), z(0) {}

  KOKKOS_INLINE_FUNCTION
  Vec3(T x_, T y_, T z_) : x(x_), y(y_), z(z_) {}

  // Vector subtraction
  KOKKOS_INLINE_FUNCTION Vec3 operator-(const Vec3& other) const {
    return Vec3(x - other.x, y - other.y, z - other.z);
  }

  // Vector addition
  KOKKOS_INLINE_FUNCTION Vec3 operator+(const Vec3& other) const {
    return Vec3(x + other.x, y + other.y, z + other.z);
  }

  // Scalar multiplication
  KOKKOS_INLINE_FUNCTION Vec3 operator*(T scalar) const {
    return Vec3(x * scalar, y * scalar, z * scalar);
  }

  // Access operator
  KOKKOS_INLINE_FUNCTION T operator[](int i) const {
    if(i == 0) return x;
    if(i == 1) return y;
    return z;
  }

  // Dot product
  KOKKOS_INLINE_FUNCTION T dot(const Vec3& other) const {
    return x * other.x + y * other.y + z * other.z;
  }

  // Euclidean norm
  KOKKOS_INLINE_FUNCTION T norm() const {
    return std::sqrt(x * x + y * y + z * z);
  }

  // Euclidean norm squared
  KOKKOS_INLINE_FUNCTION T norm2() const {
    return x * x + y * y + z * z;
  }

  // cross product
  KOKKOS_INLINE_FUNCTION Vec3 cross(const Vec3& other) const {
    return Vec3(y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x);
  }

};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Vec3<T> & v){
    os<<"Vec: ("<<v.x<<","<<v.y<<","<<v.z<<")";
    return os;
}


#endif

