#include "../include/rw_functions.hpp"
#include "../include/geom_functions.hpp"

#include <algorithm>
#include <sstream>
#include <fstream>
#include <iterator>
#include <list>
#include <tuple>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Polyhedron_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Polyhedron_3<K>                                 Polyhedron_3;
typedef K::Plane_3                                            Plane;
typedef K::Point_3                                            Point_3;



std::tuple<double, double, double> compute_plane_equation(const Polyhedron_3::Face& face) {
    auto halfedge = face.halfedge();
    Point_3 p1 = halfedge->vertex()->point();
    Point_3 p2 = halfedge->next()->vertex()->point();
    Point_3 p3 = halfedge->next()->next()->vertex()->point();

    K::Vector_3 v1 = p2 - p1;
    K::Vector_3 v2 = p3 - p1;

    K::Vector_3 normal = CGAL::cross_product(v1, v2);
    double A = normal.x();
    double B = normal.y();
    double C = normal.z();
    double length = std::sqrt(A * A + B * B + C * C);
    A /= length;
    B /= length;
    C /= length;

    return std::make_tuple(A, B, C);
}


std::tuple<int, int, int> read_spoil_structures_from_file(std::istream& in
        , std::vector<Point3D>& arr_points3d, std::vector<Point3D>& arr_norm3d
        , std::vector<Edge>& arr_edges3d) {


    std::list<Plane> planes;
    Polyhedron_3 poly;

#if 0
    std::string line;
    std::string key("#  indices of vertices");
    std::string stop("#");
    for(;;) {
        std::getline(in, line);
        if (line.find(key) != std::string::npos) 
            break;
    }
    std::getline(in, line);
    int counter = 1;
    while(line.find(stop) == std::string::npos) {
        std::getline(in, line);
        double a, b, c, d;
        sscanf(line.c_str(), "%lf %lf %lf %lf", &a, &b, &c, &d);
        if (counter == 0) {
            typename K::Plane_3 plane(a, b, c, d);
            planes.push_back(plane);
        } else {
            typename K::Plane_3 plane(a, b, c, d);
            planes.push_back(plane);
        }
        std::getline(in, line);
        ++counter;
    }
#endif
    std::string line;
    //std::string key("# plane coeff");
    std::string stop("#");
    std::string key("indices_of_vertices");
    for(;;) {
        std::getline(in, line);
        if (line.find(key) != std::string::npos) 
            break;
    }
    std::getline(in, line);
    int counter = 1;
    srand(time(NULL));
    while(line.find(stop) == std::string::npos) {
        double a, b, c, d;
        int unused1, unused2;
        sscanf(line.c_str(), "%d %d %lf %lf %lf %lf", &unused1, &unused2, &a, &b, &c, &d);
        if (c > 0) {
            typename K::Plane_3 plane(a, b, c, d + 0.5);
            planes.push_back(plane);
        } else {
            typename K::Plane_3 plane(a, b, c, d);
            planes.push_back(plane);
        }
        std::getline(in, line);
        std::getline(in, line);
        ++counter;
    }

    CGAL::halfspace_intersection_3(planes.begin(), planes.end(), poly);
    for(auto it = poly.vertices_begin(); it != poly.vertices_end(); ++it) {
        arr_points3d.emplace_back((*it).point().x(), (*it).point().y(), (*it).point().z());
    }
    for(auto it = poly.facets_begin(); it != poly.facets_end(); ++it) {
        auto [A, B, C] = compute_plane_equation(*it);
        arr_norm3d.emplace_back(A, B, C);
    }

    for(auto edge = poly.edges_begin(); edge != poly.edges_end(); ++edge) {
        auto face1 = edge->face();
        auto face2 = edge->opposite()->face();

        auto [A1, B1, C1] = compute_plane_equation(*face1);
        Point3D norm1(A1, B1, C1);
        auto [A2, B2, C2] = compute_plane_equation(*face2);
        Point3D norm2(A2, B2, C2);

        Point3D p1(edge->vertex()->point().x(), edge->vertex()->point().y(), edge->vertex()->point().z());
        Point3D p2(edge->opposite()->vertex()->point().x(), edge->opposite()->vertex()->point().y()
                , edge->opposite()->vertex()->point().z());
        int point_index_1 = 0;
        int point_index_2 = 0;
        int norm_index_1 = 0;
        int norm_index_2 = 0;
        for(int i = 0; i < arr_points3d.size(); ++i) {
            if (arr_points3d[i] == p1) {
                point_index_1 = i;
            }
            if(arr_points3d[i] == p2) {
                point_index_2 = i;
            }
        }
        for(int i = 0; i < arr_norm3d.size(); ++i) {
            if (arr_norm3d[i] == norm1) {
                norm_index_1 = i;
            }
            if(arr_norm3d[i] == norm2) {
                norm_index_2 = i;
            }
        }
        Edge ed{point_index_1, point_index_2, norm_index_1, norm_index_2};
        arr_edges3d.push_back(ed);
    }
    return std::make_tuple(arr_points3d.size(), arr_norm3d.size(), arr_edges3d.size());
}

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
#if 0
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
#endif
    for (int i = 0; i < num_facets; ++i) {
        Point3D point;
        std::getline(in, line);
        std::getline(in, line);
        std::istringstream iss(line);
        iss >> point.x >> point.y >> point.z;
        arr_norm3d.push_back(point);
    }




    for(int i = 0; i < 2; ++i) {
        std::getline(in, line);
        std::cout << line << std::endl;
    }
#if 0
    for(int i = 0; i < num_edges; ++i) {
        Edge ed;
        std::getline(in, line);
        std::istringstream iss(line);
        iss >> ed.vert1_id >> ed.vert2_id >> ed.facet1_id >> ed.facet2_id;
        arr_edges3d.push_back(ed);
    }
#endif
    for(int i = 0; i < num_edges; ++i) {
        Edge ed;
        std::getline(in, line);
        std::istringstream iss(line);
        double ignore;
        iss >> ignore;
        iss >> ed.facet1_id >> ed.facet2_id >> ed.vert1_id >> ed.vert2_id;
        arr_edges3d.push_back(ed);
    }

    return std::make_tuple(num_vertices, num_facets, num_edges);
}

void read_points_in_contour(boost::geometry::model::linestring<Point2D>& in_contour
        , std::ifstream& in) {
    std::istream_iterator<Point2D> ii(in);
    Point2D point_1 = *ii;
    in_contour.push_back(point_1);
    ++ii;

    for(Point2D point_2; ii != std::istream_iterator<Point2D>{}; ++ii) {
        point_2 = *ii;
        double dist = boost::geometry::distance(point_1, point_2);
        if (dist >= EPSILON) {
        	in_contour.push_back(point_2);
        }
        point_1 = point_2;
    }
}

