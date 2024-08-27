#pragma once 

#include "point2d.h"
#include "point3d.h"
#include "edge.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>


BOOST_GEOMETRY_REGISTER_POINT_2D(Point2D, double, boost::geometry::cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_3D(Point3D, double, boost::geometry::cs::cartesian, x, y, z);

constexpr double EPSILON = 1e-13;
inline constexpr int PRECISION = 12;
inline constexpr double BOTTOM_ERR = 1e-1;

void make_projection(double angle, std::vector<Point3D>& arr_points3d
        , boost::geometry::model::multi_point<Point2D, std::vector>& mpoint);
Point2D find_minimal_y(boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull);
std::pair<double, double> find_borders_for_x(Point2D lowest
        , boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull);
void projection_convex_hull(std::vector<Point3D>& arr_points3d
        , boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull, double angle);
void projection_square_complexity(std::vector<Point3D>& arr_points3d
        , std::vector<Point3D>& arr_norm3d
        , std::vector<Edge>& arr_edges3d
        , boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull, double angle); 




