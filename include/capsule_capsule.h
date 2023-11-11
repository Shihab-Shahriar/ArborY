#include "segment_segment.h"
#include "vector.h"

/**
 * @brief Compute the distance between two capsules
 * 
 * @tparam Scalar float or double
 * @tparam Vec3 3D vector type
 * 
 * @param capsule0_l1 First endpoint of the first capsule
 * @param capsule0_l2 Second endpoint of the first capsule
 * @param capsule0_radius Radius of the first capsule
 * @param capsule1_l1 First endpoint of the second capsule
 * @param capsule1_l2 Second endpoint of the second capsule
 * @param capsule1_radius Radius of the second capsule
 * 
 * @return Scalar Distance between the two capsules, 0 if overlapping
*/
template <typename Scalar, typename Vec3=Vec3<Scalar>>
Scalar CapsuleCapsuleDistance(Vec3 capsule0_l1, Vec3 capsule0_l2, double capsule0_radius, 
                            Vec3 capsule1_l1, Vec3 capsule1_l2, double capsule1_radius) {
    auto distance = SegmentSegmentDistanceSquared<Scalar>(capsule0_l1, capsule0_l2, capsule1_l1, capsule1_l2);
    double radiusSum = capsule0_radius + capsule1_radius;
    return Kokkos::fmax(Kokkos::sqrt(distance) - radiusSum, 0.0);
}
