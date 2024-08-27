#include "../include/functions.h"

#include <algorithm>
#include <sstream>
#include <iterator>

std::tuple<int, int, int> read_structures_from_file(std::istream& in
        , std::vector<Point3D>& arr_points3d, std::vector<Point3D>& arr_norm3d
        , std::vector<Edge>& arr_edges3d) {
    int num_vertices, num_facets, num_edges;

    std::string line;
    for(int i = 0; i < 3; ++i) { 
        std::getline(in, line);
        std::cout << line << std::endl;
    }

	std::istringstream iss(line);
    iss >> num_vertices;
    iss >> num_facets;
    iss >> num_edges;

    for(int i = 0; i < 2; ++i) {
        std::getline(in, line);
        std::cout << line << std::endl;
    }

    std::copy_n(std::istream_iterator<Point3D>(in), num_vertices, std::back_inserter(arr_points3d));

    for(int i = 0; i < 4; ++i) {
        std::getline(in, line);
        std::cout << line << std::endl;
    }

    for (int i = 0; i < num_facets; ++i) {
        Point3D point;
        std::getline(in, line);
        std::istringstream iss(line);
        double ignore;
        iss >> ignore >> ignore;
        iss >> point.x >> point.y >> point.z;
        arr_norm3d.push_back(point);
        std::getline(in, line);
    }

    for(int i = 0; i < 2; ++i) {
        std::getline(in, line);
        std::cout << line << std::endl;
    }

    for(int i = 0; i < num_edges; ++i) {
        Edge ed;
        std::getline(in, line);
        std::istringstream iss(line);
        iss >> ed.vert1_id >> ed.vert2_id >> ed.facet1_id >> ed.facet2_id;
        arr_edges3d.push_back(ed);
    }

    return std::make_tuple(num_vertices, num_facets, num_edges);
}

