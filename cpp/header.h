#pragma once

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <ranges>
#include <algorithm>

#include <random>

#include <boost/geometry.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/math/constants/constants.hpp>

#include <boost/geometry/geometries/register/point.hpp>

#include <boost/geometry/geometries/geometries.hpp>


inline constexpr int POINTS_BETWEEN = 30;

struct Point {
    double x = 0, y = 0;
    Point() = default;
    Point(double x, double y) noexcept : x(x), y(y) {}
    Point(const Point& p) noexcept : x(p.x), y(p.y) {}
    Point(Point&& p): x(p.x), y(p.y) {}
    Point& operator = (const Point& p) noexcept {
        x = p.x; y = p.y;
        return *this;
    }
    Point& operator = (Point&& p) noexcept {
        x = p.x; y = p.y;
        return *this;
    }
    friend Point operator + (const Point& p1, const Point& p2) noexcept {
        return Point(p1.x + p2.x, p1.y + p2.y);
    }
    friend Point operator - (const Point& p1, const Point& p2) noexcept {
        return Point(p1.x - p2.x, p1.y - p2.y);
    }
    friend Point operator / (const Point& p, double del) {
        return Point(p.x / del, p.y / del);
    }
    friend Point operator * (double mul, const Point& p) {
        return Point(p.x * mul, p.y * mul);
    }
    friend Point operator * (const Point& p, double mul) {
        return Point(p.x * mul, p.y * mul);
    }
    ~Point() = default;
};

std::istream& operator >> (std::istream& is, Point& p);
std::ostream& operator << (std::ostream& os, const Point& p);


BOOST_GEOMETRY_REGISTER_POINT_2D(Point, double, boost::geometry::cs::cartesian, x, y)

Point rotate(Point p, double angle);
void grind(std::vector<Point>& v, Point point_1, Point point_2);
void read_points_in_vector(std::vector<Point>& v, std::ifstream& in);
void spoil_and_get_protrusions(std::vector<Point>& v, std::vector<boost::geometry::model::polygon<Point, false, true, std::vector>>& vec_of_pol); 
std::pair<double, double> write_intersection(std::vector<boost::geometry::model::polygon<Point>>& dif
        , std::ofstream& out_intersection); 



