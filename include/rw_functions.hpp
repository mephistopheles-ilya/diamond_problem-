#pragma once

#include "point3d.hpp"
#include "point2d.hpp"
#include "edge.hpp"

#include <tuple>
#include <vector>

std::tuple<int, int, int> read_structures_from_file(std::istream& in
        , std::vector<Point3D>& arr_points3d, std::vector<Point3D>& arr_norm3d
        , std::vector<Edge>& arr_edges3d);
void read_points_in_contour(boost::geometry::model::linestring<Point2D>& in_contour, std::ifstream& in);
std::tuple<int, int, int> read_spoil_structures_from_file(std::istream& in
        , std::vector<Point3D>& arr_points3d, std::vector<Point3D>& arr_norm3d
        , std::vector<Edge>& arr_edges3d, bool is_change_shift,  double parametr_of_shift
        , bool is_change_azimuth, double parametr_of_azimuth, bool is_change_slope, double parametr_of_slope
        , int sign_of_c, std::string distribution, double sigma);

