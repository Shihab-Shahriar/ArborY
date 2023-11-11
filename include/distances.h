#pragma once

#include <Kokkos_Core.hpp>
#include "vector.h"
#include <iostream>

/**
 * 
 */

using Vec3d = Vec3<double>;

// A simple point-to-point distance function as an example.
KOKKOS_INLINE_FUNCTION
double PointToPointDistance(const Vec3d& p1, const Vec3d& p2) {
  return Kokkos::sqrt((p2.x - p1.x) * (p2.x - p1.x) +
              (p2.y - p1.y) * (p2.y - p1.y) +
              (p2.z - p1.z) * (p2.z - p1.z));
}

// Point to line segment distance
KOKKOS_INLINE_FUNCTION
double PointToLineSegmentDistance(const Vec3d& p, const Vec3d& l1, const Vec3d& l2) {
  Vec3d v = l2 - l1;
  Vec3d w = p - l1;

  double c1 = w.dot(v);
  double c2 = v.dot(v);
  double t = c1 / c2;

  Vec3d distanceVec;
  if (t < 0) {
    distanceVec = p - l1;
  } else if (t > 1) {
    distanceVec = p - l2;
  } else {
    Vec3d projection = l1 + v * t;
    distanceVec = p - projection;
  }

  return distanceVec.norm();
}      

// Point to plane distance
KOKKOS_INLINE_FUNCTION
double PointToPlaneDistance(const Vec3d& p, const Vec3d& planeNormal, double planeConstant) {
  double distance = p.dot(planeNormal) + planeConstant;
  return distance / planeNormal.norm();
}

// Point to Line distance
KOKKOS_INLINE_FUNCTION
double PointToLineDistance(const Vec3d& p, const Vec3d& l1, const Vec3d& l2) {
  Vec3d v = l2 - l1;
  Vec3d w = p - l1;

  double c1 = w.dot(v);
  double c2 = v.dot(v);
  double t = c1 / c2;

  Vec3d distanceVec;
  if (t < 0) {
    distanceVec = p - l1;
  } else if (t > 1) {
    distanceVec = p - l2;
  } else {
    Vec3d projection = l1 + v * t;
    distanceVec = p - projection;
  }

  return distanceVec.norm();
}

// point to sphere distance
KOKKOS_INLINE_FUNCTION
double PointToSphereDistance(const Vec3d& p, const Vec3d& center, double radius) {
  double distance = PointToPointDistance(p, center);
  return Kokkos::fmax(distance - radius, 0.0);
}

// point to box distance
KOKKOS_INLINE_FUNCTION
double PointToBoxDistance(const Vec3d& p, const Vec3d& min, const Vec3d& max) {
  using T = double;
  T distanceVec[3];
  for (int i = 0; i < 3; i++) {
    if (p[i] < min[i]) {
      distanceVec[i] = min[i] - p[i];
    } else if (p[i] > max[i]) {
      distanceVec[i] = p[i] - max[i];
    } else {
      distanceVec[i] = 0.0;
    }
  }
  Vec3d d(distanceVec[0], distanceVec[1], distanceVec[2]);
  return d.norm();
}


/*
 * Copied code from GeometricTools: 
 * https://github.com/davideberly/GeometricTools/blob/master/GTE/Mathematics/DistPointTriangle.h
 * 
 * Explanation of what's going on:
 * https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
 * 
 * Boost Software License 1.0 License
 */
KOKKOS_INLINE_FUNCTION
double PointToTriangleDistance(const Vec3d& p, const Vec3d& t0, const Vec3d& t1, const Vec3d& t2) {
  Vec3d diff = t0 - p;
  Vec3d edge0 = t1 - t0;
  Vec3d edge1 = t2 - t0;

  using T = double;

  T zero = 0.0;
  T one = 1.0;
  T two = 2.0;

  T a00 = edge0.dot(edge0);
  T a01 = edge0.dot(edge1);
  T a11 = edge1.dot(edge1);
  T b0 = diff.dot(edge0);
  T b1 = diff.dot(edge1);
  T det = Kokkos::fmax(a00 * a11 - a01 * a01, 0.0);
  T s = a01 * b1 - a11 * b0;
  T t = a01 * b0 - a00 * b1;

  if (s + t <= det)
  {
      if (s < 0.0)
      {
          if (t < 0.0)  // region 4
          {
              if (b0 < 0.0)
              {
                  t = 0.0;
                  if (-b0 >= a00)
                  {
                      s = 1.0;
                  }
                  else
                  {
                      s = -b0 / a00;
                  }
              }
              else
              {
                  s = 0.0;
                  if (b1 >= 0.0)
                  {
                      t = 0.0;
                  }
                  else if (-b1 >= a11)
                  {
                      t = 1.0;
                  }
                  else
                  {
                      t = -b1 / a11;
                  }
              }
          }
          else  // region 3
          {
              s = zero;
              if (b1 >= zero)
              {
                  t = zero;
              }
              else if (-b1 >= a11)
              {
                  t = one;
              }
              else
              {
                  t = -b1 / a11;
              }
          }
      }
      else if (t < zero)  // region 5
      {
          t = zero;
          if (b0 >= zero)
          {
              s = zero;
          }
          else if (-b0 >= a00)
          {
              s = one;
          }
          else
          {
              s = -b0 / a00;
          }
      }
      else  // region 0
      {
          // minimum at interior point
          s /= det;
          t /= det;
      }
  }
  else
  {
      T tmp0{}, tmp1{}, numer{}, denom{};

      if (s < zero)  // region 2
      {
          tmp0 = a01 + b0;
          tmp1 = a11 + b1;
          if (tmp1 > tmp0)
          {
              numer = tmp1 - tmp0;
              denom = a00 - two * a01 + a11;
              if (numer >= denom)
              {
                  s = one;
                  t = zero;
              }
              else
              {
                  s = numer / denom;
                  t = one - s;
              }
          }
          else
          {
              s = zero;
              if (tmp1 <= zero)
              {
                  t = one;
              }
              else if (b1 >= zero)
              {
                  t = zero;
              }
              else
              {
                  t = -b1 / a11;
              }
          }
      }
      else if (t < zero)  // region 6
      {
          tmp0 = a01 + b1;
          tmp1 = a00 + b0;
          if (tmp1 > tmp0)
          {
              numer = tmp1 - tmp0;
              denom = a00 - two * a01 + a11;
              if (numer >= denom)
              {
                  t = one;
                  s = zero;
              }
              else
              {
                  t = numer / denom;
                  s = one - t;
              }
          }
          else
          {
              t = zero;
              if (tmp1 <= zero)
              {
                  s = one;
              }
              else if (b0 >= zero)
              {
                  s = zero;
              }
              else
              {
                  s = -b0 / a00;
              }
          }
      }
      else  // region 1
      {
          numer = a11 + b1 - a01 - b0;
          if (numer <= zero)
          {
              s = zero;
              t = one;
          }
          else
          {
              denom = a00 - two * a01 + a11;
              if (numer >= denom)
              {
                  s = one;
                  t = zero;
              }
              else
              {
                  s = numer / denom;
                  t = one - s;
              }
          }
      }
  }

  Vec3d projected = t0 + edge0 * s + edge1 * t;
  diff = p - projected;
  return diff.norm();
}


/*
 * @brief: Point to capsule distance (also known as spherocylinder)
 * 
 * @param p: point
 * @param l1: Cylinder axis start 
 * @param l2: Cylinder axis end
 * @param radius: radius of cylinder
 * 
 * @return: distance from point to capsule
 * 
 */
KOKKOS_INLINE_FUNCTION
double PointToCapsuleDistance(const Vec3d& p, const Vec3d& l1, const Vec3d& l2, double radius) {
  Vec3d axis = l2 - l1;
  Vec3d w = p - l1;

  double c1 = w.dot(axis);
  double c2 = axis.dot(axis);
  double t = c1 / c2;

  //printf("t: %f\n", t);

  // Algo taken from here: https://liris.cnrs.fr/Documents/Liris-1297.pdf
  Vec3d h = l1 + axis * t;
  Vec3d y = p - h;

  if(t>=0.0 && t<=1.0){ // projection is on line segment
    return Kokkos::fmax(y.norm() - radius, 0.0);
  }

  
  // projection is outside line segment
  return Kokkos::fmin(PointToSphereDistance(p, l1, radius), PointToSphereDistance(p, l2, radius));
}


/**
 * @brief: AABB to AABB distance
 * 
 * @param min1: min corner of AABB 1
 * @param max1: max corner of AABB 1
 * @param min2: min corner of AABB 2
 * @param max2: max corner of AABB 2
 * 
 * @return: distance between AABBs
*/
KOKKOS_INLINE_FUNCTION
double BoxToBoxDistance(const Vec3d& min1, const Vec3d& max1, const Vec3d& min2, const Vec3d& max2) {
  using T = double;
  T distanceVec[3];
  for (int i = 0; i < 3; i++) {
    if (min1[i] > max2[i]) {
      distanceVec[i] = min1[i] - max2[i];
    } else if (min2[i] > max1[i]) {
      distanceVec[i] = min2[i] - max1[i];
    } else {
      distanceVec[i] = 0.0;
    }
  }
  Vec3d d(distanceVec[0], distanceVec[1], distanceVec[2]);
  return d.norm();
}

/**
 * @brief: AABB to sphere distance
 * 
 * @param min: min corner of AABB
 * @param max: max corner of AABB
 * @param center: center of sphere
 * @param radius: radius of sphere
 * 
 * @return: distance between AABB and sphere
*/
KOKKOS_INLINE_FUNCTION
double BoxToSphereDistance(const Vec3d& min, const Vec3d& max, const Vec3d& center, double radius) {
  double distance = PointToBoxDistance(center, min, max);
  return Kokkos::fmax(distance - radius, 0.0);
}

