#pragma once

#include "point3d.h"
#include "edge.h"

#include <tuple>
#include <vector>

std::tuple<int, int, int> read_structures_from_file(std::istream& in
        , std::vector<Point3D>& arr_points3d, std::vector<Point3D>& arr_norm3d
        , std::vector<Edge>& arr_edges3d);

