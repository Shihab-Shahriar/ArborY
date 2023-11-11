#include <CGAL/Simple_cartesian.h>
#include <CGAL/Plane_3.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 Plane_3;

void pointToPoint(){
    Point_3 p1(1.0, 2.0, 3.0);
    Point_3 p2(4.0, 5.0, 6.0);

    double distance = CGAL::squared_distance(p1, p2);

    std::cout << "The sq distance from the point to the point is: " << distance << std::endl;
}

void pointToSegment(){
    typedef Kernel::Segment_3 Segment_3;

    Point_3 point(1.0, 2.0, 3.0); 

    Point_3 seg_point_1(0.0, -1.0, 10.0); 
    Point_3 seg_point_2(1.0, 3.0, -9.0); // Second endpoint 
    Segment_3 segment(seg_point_1, seg_point_2);

    double distance = CGAL::squared_distance(point, segment);

    std::cout << "The distance from the point to the segment is: " << distance << std::endl;
}

void pointToPlane(){
    Point_3 point(1.0, 2.0, 3.0); 
    Plane_3 plane(1.0, 2.0, 3.0, 4.0);

    double distance = CGAL::squared_distance(point, plane);

    std::cout << "The distance from the point to the plane is: " << distance << std::endl;

}

void pointToTriangle(){
    typedef Kernel::Triangle_3 Triangle_3;

    Point_3 point(1.0, 2.0, 3.0); 

    Point_3 tri_point_1(0.0, -1.0, 2.0); 
    Point_3 tri_point_2(1.0, 3.0, 2.0); 
    Point_3 tri_point_3(2.0, 9.0, 2.0); 
    Triangle_3 triangle(tri_point_1, tri_point_2, tri_point_3);

    double distance = CGAL::squared_distance(point, triangle);

    std::cout << "The distance from the point to the triangle is: " << distance << std::endl;
}

void pointToLine(){
    typedef Kernel::Line_3 Line_3;

    Point_3 point(1.0, 2.0, 3.0); 

    Point_3 line_point_1(0.0, -1.0, 10.0); 
    Point_3 line_point_2(1.0, 3.0, -9.0); // Second endpoint 
    Line_3 line(line_point_1, line_point_2);

    double distance = CGAL::squared_distance(point, line);

    std::cout << "The distance from the point to the line is: " << distance << std::endl;
}

void segmentToSegment(){
    typedef Kernel::Segment_3 Segment_3;

    Point_3 seg1_point_1(0.0, 0.0, 0.0); 
    Point_3 seg1_point_2(2.0, 0.0, 0.0); // Second endpoint 
    Segment_3 segment1(seg1_point_1, seg1_point_2);

    Point_3 seg2_point_1(0.0, 1.0, 0.0); 
    Point_3 seg2_point_2(2.0, 1.0, 0.0); // Second endpoint 
    Segment_3 segment2(seg2_point_1, seg2_point_2);

    Point_3 seg3_point_1(7.0, 11.0, 9.0);
    Point_3 seg3_point_2(3.0, -8.0, 12.0); 
    Segment_3 segment3(seg3_point_1, seg3_point_2);

    double distance = CGAL::squared_distance(segment1, segment2);
    double distance2 = CGAL::squared_distance(segment1, segment3);

    std::cout << "The distance from the segment to the segment is: " << distance << std::endl;
    std::cout << "The distance from the segment to the segment is: " << distance2 << std::endl;
}

int main() {
    pointToPoint();
    pointToSegment();
    pointToPlane();
    pointToTriangle();
    pointToLine();
    segmentToSegment();

    return 0;
}
