#include "distances.h"
#include "segment_segment.h"
#include "capsule_capsule.h"
#include "vector.h"

#include <Kokkos_Core.hpp>

using Vec3d = Vec3<double>;

void point_point(){
    Vec3d p1(1.0, 2.0, 3.0);
    Vec3d p2(4.0, 5.0, 6.0);

    double distance = PointToPointDistance(p1, p2);

    std::cout << "The distance from the point to point is: " << distance*distance << std::endl;
}

void point_segment(){
    Vec3d point(1.0, 2.0, 3.0); 

    Vec3d seg_point_1(0.0, -1.0, 10.0); 
    Vec3d seg_point_2(1.0, 3.0, -9.0); // Second endpoint 
    
    double distance = PointToLineSegmentDistance(point, seg_point_1, seg_point_2);

    std::cout << "The distance from the point to the segment is: " << distance*distance << std::endl;
}

void point_line(){
    Vec3d point(1.0, 2.0, 3.0); 

    // the line passes through these two points
    Vec3d point_1(0.0, -1.0, 10.0); 
    Vec3d point_2(1.0, 3.0, -9.0); // Second endpoint 
    
    double distance = PointToLineDistance(point, point_1, point_2);

    std::cout << "The distance from the point to the line is: " << distance*distance << std::endl;
}

void point_plane(){
    Vec3d point(1.0, 2.0, 3.0); 
    Vec3d normal(1.0, 2.0, 3.0);
    double constant = 4.0;

    double distance = PointToPlaneDistance(point, normal, constant);

    std::cout << "The distance from the point to the plane is: " << distance*distance << std::endl;
}

void point_triangle(){
    //TODO: test when point inside triangle

    Vec3d point(1.0, 2.0, 3.0); 

    Vec3d tri_point_1(0.0, -1.0, 2.0); 
    Vec3d tri_point_2(1.0, 3.0, 2.0); 
    Vec3d tri_point_3(2.0, 9.0, 2.0); 

    double distance = PointToTriangleDistance(point, tri_point_1, tri_point_2, tri_point_3);

    std::cout << "The distance from the point to the triangle is: " << distance*distance << std::endl;
}


void point_capsule(){ //spherocylinder
    Vec3d l1(0.0, 0.0, 0.0);
    Vec3d l2(2.0, 0.0, 0.0);
    double radius = 0.5;

    Vec3d point(3.0, 0.0, 0.0);

    double distance = PointToCapsuleDistance(point, l1, l2, radius);
    std::cout<<"The distance from the point to the capsule is: "<<distance<<std::endl;

    point = Vec3d(1.0, 0.0, 0.0);
    distance = PointToCapsuleDistance(point, l1, l2, radius);
    std::cout<<"The distance from the point to the capsule is: "<<distance<<std::endl;

    point = Vec3d(1.0, 2.0, 0.0);
    distance = PointToCapsuleDistance(point, l1, l2, radius);
    std::cout<<"The distance from the point to the capsule is: "<<distance<<std::endl;
}

void segment_segment(){
    Vec3d l1_1(0.0, 0.0, 0.0);
    Vec3d l1_2(2.0, 0.0, 0.0);

    Vec3d l2_1(0.0, 1.0, 0.0);
    Vec3d l2_2(2.0, 1.0, 0.0);


    Vec3d l3_1(7.0, 11.0, 9.0);
    Vec3d l3_2(3.0, -8.0, 12.0); 

    double distance = SegmentSegmentDistanceSquared<double, Vec3d>(l1_1, l1_2, l2_1, l2_2);
    double distance2 = SegmentSegmentDistanceSquared<double, Vec3d>(l1_1, l1_2, l3_1, l3_2);


    std::cout<<"The sq distance from the segment to the segment is: "<<distance<<std::endl;
    std::cout<<"The sq distance from the segment to the segment is: "<<distance2<<std::endl;
}

void capsule_capsule(){
    Vec3d l1_1(0.0, 0.0, 0.0);
    Vec3d l1_2(2.0, 0.0, 0.0);
    double radius1 = 0.5;

    Vec3d l2_1(0.0, 1.0, 0.0);
    Vec3d l2_2(2.0, 1.0, 0.0);
    double radius2 = 0.45;

    Vec3d l3_1(7.0, 11.0, 9.0);
    Vec3d l3_2(3.0, -8.0, 12.0); 
    double radius3 = 0.75;

    // These two capsules are parallel, same plane, easy to visualize 
    double distance = CapsuleCapsuleDistance<double, Vec3d>(l1_1, l1_2, radius1, l2_1, l2_2, radius2);
    
    double distance2 = CapsuleCapsuleDistance<double, Vec3d>(l1_1, l1_2, radius1, l3_1, l3_2, radius3);

    std::cout<<"The distance from the capsule to the capsule is: "<<distance<<std::endl;
    std::cout<<"The distance from the capsule to the capsule is: "<<distance2<<std::endl;

}


int main(int argc, char *argv[]){
    Kokkos::ScopeGuard guard(argc, argv);

    point_point();
    point_segment();
    point_plane();
    point_triangle();
    point_line();
    point_capsule();
    segment_segment();
    capsule_capsule();

    return 0;
}