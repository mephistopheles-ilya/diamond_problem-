#include <iostream>
#include <fstream>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>

#include "include/point2d.hpp"
#include "include/segment2d.hpp"
#include "include/geom_functions.hpp"

using line_string = boost::geometry::model::linestring<Point2D>;
using segment = boost::geometry::model::segment<Point2D>;
using multy_point = boost::geometry::model::multi_point<Point2D>;

inline constexpr double x_left = -100;
inline constexpr double x_right = 100;
inline constexpr double max = 200;

#if 0
double compare_contours(line_string& points, line_string& polyline, std::vector<Segment2D>& segs) {
    double y_max = points[0].y;
    double y_min = std::max(points[0].y, points[points.size() - 1].y);
    for(auto& p: polyline) {
        if(p.y > y_max) y_max = p.y;
    }
    double error = 0;
    double delta_x = 0;
    int index_1 = 0;
    int index_2 = 1;
    Point2D point1 = polyline[0];
    Point2D point2 = polyline[1];
    double y = points[index_1].y;

    while(y < point1.y) {
        ++index_1;
        if (index_1 < points.size()) return error;
        y = points[index_1].y;
    }

    for(;;) {
        double y_l = point1.y;
        double y_h = point2.y;
        if(y_l > y_h) std::swap(y_l, y_h);
        while(!(y >= y_l && y <= y_h)) {
            ++index_2;
            if(index_2 > polyline.size()) return error;
            point1 = point2;
            point2 = polyline[index_2];
            y_l = point1.y;
            y_h = point2.y;
            if(y_l > y_h) std::swap(y_l, y_h);
        }

        Point2D vec = point2 - point1;
        double A = -vec.y;
        double B = vec.x;
        double C = point1.x * vec.y - point1.y * vec.x;
        Point2D intersection(-(B * y + C) / A, y); 
        delta_x = std::fabs(points[index_1].x - intersection.x);

        Segment2D seg{intersection, points[index_1]};
        segs.push_back(seg);

        error += delta_x;
        ++index_1;
        if (index_1 > points.size()) return error;
        y = points[index_1].y;
        while(y > y_max) {
            ++index_1;
            if(index_1 > points.size()) return error;
            y = points[index_1].y;
        }
        while(y < y_min) {
            ++index_1;
            if(index_1 > points.size()) return error;
            y = points[index_1].y;
        }
    }
    return error;
}
#endif


double compare_contours(line_string& points, line_string& polyline, std::vector<Segment2D>& segs) {
    double y_max = points[0].y;
    double y_min = std::max(points[0].y, points[points.size() - 1].y);
    for(auto& p: polyline) {
        if(p.y > y_max) y_max = p.y;
    }

    line_string line;
    double error = 0;
    std::vector<Point2D> intersections;
    for(auto& p: points) {
        if(p.y <= y_max && p.y >= y_min) {
            Point2D left(x_left, p.y);
            Point2D right(x_right, p.y);
            line.push_back(left);
            line.push_back(right);
            if(boost::geometry::intersection(polyline, line, intersections) == true) {
                double delta_x = max;
                Point2D intersection;
                for(auto p_int: intersections) {
                    if(std::fabs(p_int.x - p.x) < delta_x) {
                        delta_x = std::fabs(p_int.x - p.x);
                        intersection = p_int;
                    }
                }
                error += delta_x;
                Segment2D seg{p, intersection};
                segs.push_back(seg);
                intersections.clear();
            }
            line.clear();
        }
    }
    return error;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage <directory1> <directory2> <amount_of_files>\n";
        return -1;
    }
    int amount_of_files = 0;
    if (sscanf(argv[3], "%d", &amount_of_files) != 1) {
        std::cout << "Wrong amount of files" << std::endl;
        return -1;
    }
    for(int i = 0; i <= amount_of_files; ++i) {
        std::string contour1 = std::string(argv[1]) + std::string("/contour_") + std::to_string(i);
        std::string contour2 = std::string(argv[2]) + std::string("/contour_") + std::to_string(i);
        std::ifstream in1(contour1);
        std::ifstream in2(contour2);
        std::istream_iterator<Point2D> it1(in1);
        std::istream_iterator<Point2D> it2(in2);

        line_string line1(it1, std::istream_iterator<Point2D>{});
        line_string line2(it2, std::istream_iterator<Point2D>{});
        assert(("No points in file1", line1.size() > 3));
        assert(("No points in file2", line2.size() > 3));

        std::vector<Segment2D> segs;
        double error = compare_contours(line1, line2, segs);
        std::cout << error << std::endl;
    }
    return 0;
}

