#pragma once 

#include "point2d.h"
#include "point3d.h"
#include "edge.h"
#include "segment2d.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>

#include <unordered_map>
#include <unordered_set>


BOOST_GEOMETRY_REGISTER_POINT_2D(Point2D, double, boost::geometry::cs::cartesian, x, y);
BOOST_GEOMETRY_REGISTER_POINT_3D(Point3D, double, boost::geometry::cs::cartesian, x, y, z);

constexpr double EPSILON = 1e-13;
inline constexpr int PRECISION = 12;
inline constexpr double BOTTOM_ERR = 1e-1;
inline constexpr double __DISTANCE_BETWEEN_POINTS = 0.008;


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




void hash_combine(size_t& seed, size_t hash_value);

struct hash_t {
    std::size_t operator()(const Point2D& p) const {
        std::size_t seed = std::hash<decltype(p.x)>{}(p.x);
        hash_combine(seed, std::hash<decltype(p.y)>{}(p.y));
        return seed;
    }
};

using graph_t = std::unordered_map<Point2D, std::vector<Point2D>, hash_t>;
using set_t = std::unordered_set<Point2D, hash_t>;

void restore_polygon(const graph_t& graph, const Point2D& point, set_t& visited,
                     std::vector<Point2D>& polygon);

void restore_polygons(const std::vector<Segment2D>& segments,
                      std::vector<std::vector<Point2D>>& polygons);

void projection_graph(std::vector<Point3D>& arr_points3d
        , std::vector<Point3D>& arr_norm3d
        , std::vector<Edge>& arr_edges3d
        , boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull, double angle);

Point2D rotate(Point2D p, double angle);
void spoil_and_get_protrusions(std::vector<Point2D>& v, std::vector<boost::geometry::model::polygon<Point2D, false, true, std::vector>>& vec_of_pol, std::vector<bool>& mask); 
void create_small_shifts(boost::geometry::model::linestring<Point2D>& pol, double sz_shift);
void uniform_grid_intersection(boost::geometry::model::linestring<Point2D>& line
        , boost::geometry::model::linestring<Point2D>& new_line
        , double dist_points = __DISTANCE_BETWEEN_POINTS);

