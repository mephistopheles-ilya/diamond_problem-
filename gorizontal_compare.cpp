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


double compare_contours(line_string& points, line_string& polyline, std::vector<Segment2D>& segs, double& max_dist, double& z_of_max_dist
        , double& sum_squared_dists, unsigned& num_dists, std::string cur_name, std::string& max_name){
    double y_max = points[0].y;
    double y_min = std::max(points[0].y, points[points.size() - 1].y);
    for(auto& p: polyline) {
        if(p.y > y_max) y_max = p.y;
        if(p.y < y_min) y_min = p.y;
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
                if (intersections.size() >= 2) {
                    double delta_x = max;
                    Point2D intersection;
                    for(auto p_int: intersections) {
                        if(std::fabs(p_int.x - p.x) < delta_x) {
                            delta_x = std::fabs(p_int.x - p.x);
                            intersection = p_int;
                        }
                    }
                    if (delta_x > max_dist) {
                        max_dist = delta_x;
                        z_of_max_dist = p.y;
                        max_name = cur_name;
                    }
                    sum_squared_dists += delta_x * delta_x;
                    num_dists += 1;

                    error += delta_x;
                    Segment2D seg{p, intersection};
                    segs.push_back(seg);
                    intersections.clear();
                }
            }
            line.clear();
        }
    }
    return error/points.size();
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
    double max_dist = 0, z_of_max_dist = 0, sum_squared_dists = 0;
    unsigned num_dists = 0;
    std::string max_dist_con_name;
    for(int i = 0; i < amount_of_files; ++i) {
        std::string angle_in_degrees = std::to_string(i * 180./amount_of_files);
        angle_in_degrees = angle_in_degrees.substr(0, angle_in_degrees.find(".") + 4);
        auto pos = angle_in_degrees.find(".");
        std::string zeros;
        while(pos < 3) {
            pos++;
            zeros += "0";
        }
        angle_in_degrees = zeros + angle_in_degrees;
        std::string contour1 = std::string(argv[1]) + std::string("/Contour") + angle_in_degrees + std::string(".txt");
        std::string contour2 = std::string(argv[2]) + std::string("/Contour") + angle_in_degrees + std::string(".txt");

        std::ifstream in1(contour1);
        std::ifstream in2(contour2);
        std::istream_iterator<Point2D> it1(in1);
        std::istream_iterator<Point2D> it2(in2);

        line_string line1(it1, std::istream_iterator<Point2D>{});
        line_string line2(it2, std::istream_iterator<Point2D>{});
        assert(("No points in file1", line1.size() > 3));
        assert(("No points in file2", line2.size() > 3));

        std::vector<Segment2D> segs;
        compare_contours(line1, line2, segs, max_dist, z_of_max_dist, sum_squared_dists, num_dists, contour1, max_dist_con_name);
        segs.clear();
        compare_contours(line2, line1, segs, max_dist, z_of_max_dist, sum_squared_dists, num_dists, contour1, max_dist_con_name);
 
    }
    double mediana = std::sqrt(sum_squared_dists / num_dists);
    std::cout << "Max in " << max_dist_con_name << std::endl;
    std::cout << "Max = " << max_dist << std::endl;
    std::cout << "z_max = " << z_of_max_dist << std::endl;
    std::cout << "mediana = " << mediana << std::endl;
    return 0;
}
