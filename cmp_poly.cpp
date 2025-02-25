#include <iostream>
#include <vector>
#include <iterator>

#include <boost/program_options.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/centroid.h>
#include <CGAL/squared_distance_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;
typedef K::Point_3                                          Point_3;
typedef K::Plane_3                                          Plane;

struct to_point {
    Point_3 operator()(std::iterator_traits<const Polyhedron::Halfedge_around_facet_circulator>::reference it) const {
        return it.vertex()->point();
    } 
};

struct norm_center {
    Point_3 norm;
    Point_3 center;
    typename Polyhedron::Face_iterator it1;
    typename Polyhedron::Face_iterator it2;
    double diff = -1;

    double operator()(const norm_center& nc) const {
        return CGAL::squared_distance(nc.norm, norm) + CGAL::squared_distance(nc.center, center);
    }
};

void write_facet_ply(Polyhedron::Facet_iterator fit, const std::string& filename) {
    std::ofstream out(filename);
    out << "ply" << std::endl;
    out << "format ascii 1.0" << std::endl;
    out << "element vertex " << fit->size() << std::endl; 
    out << "property double x" << std::endl;
    out << "property double y" << std::endl;
    out << "property double z" << std::endl;
    out << "element face 1" << std::endl;
    out << "property list uchar int vertex_indices" << std::endl;
    out << "end_header" << std::endl;

    // Write vertices
    auto pit = fit->facet_begin();
    for(size_t i = 0; i < fit->size(); ++i, ++pit) {
        out << pit->vertex()->point() << std::endl;
    }
    out << fit->size() << " ";
    for(size_t i = 0; i < fit->size(); ++i) {
        out << i << " ";
    }

    out << std::endl;

    out.close();
}

int main(int argc, char* argv[]) {
    std::string type_of_file;
    std::string file1;
    std::string file2;
    boost::program_options::options_description desc("All options");
    desc.add_options()
        ("type_of_file", boost::program_options::value<std::string>(&type_of_file)
         -> default_value("ply"), "ply or obj")
        ("file1", boost::program_options::value<std::string>(&file1)
         -> default_value("poly1.ply"), "first polyhedron")
        ("file2", boost::program_options::value<std::string>(&file2)
         -> default_value("poly1.ply"), "second polyhedron")
        ("help", "produce help message")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    Polyhedron poly1;
    Polyhedron poly2;
    if (type_of_file == "ply") {
        CGAL::IO::read_PLY(file1, poly1);
        CGAL::IO::read_PLY(file2, poly2);
    } else {
        CGAL::IO::read_OBJ(file1, poly1);
        CGAL::IO::read_OBJ(file2, poly2);
    }
    std::vector<norm_center> mer1;
    for(auto fit = poly1.facets_begin(); fit != poly1.facets_end(); ++fit) {
        auto begin = boost::make_transform_iterator(fit->facet_begin(), to_point());
        auto end = boost::make_transform_iterator(fit->facet_begin(), to_point());
        std::advance(end, fit->size() - 1);
    
        Point_3 center = CGAL::centroid(begin, end);
        auto h = (*fit).halfedge();
        Plane plane  = Plane(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
        double norm_length = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
        Point_3 norm = Point_3(plane.a() / norm_length, plane.b() / norm_length, plane.c() / norm_length);
        mer1.emplace_back(center, norm, fit);
    }
    std::vector<norm_center> mer2;
    for(auto fit = poly2.facets_begin(); fit != poly2.facets_end(); ++fit) {
        auto begin = boost::make_transform_iterator(fit->facet_begin(), to_point());
        auto end = boost::make_transform_iterator(fit->facet_begin(), to_point());
        std::advance(end, fit->size() - 1);
    
        Point_3 center = CGAL::centroid(begin, end);
        auto h = (*fit).halfedge();
        Plane plane  = Plane(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
        double norm_length = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
        Point_3 norm = Point_3(plane.a() / norm_length, plane.b() / norm_length, plane.c() / norm_length);
        mer2.emplace_back(center, norm, fit);
    }

    for(auto& el1: mer1) {
        double min_diff = std::numeric_limits<double>::max();
        double diff = 0;
        for(auto& el2: mer2) {
            diff = el1(el2);
            if (diff < min_diff) {
                min_diff = diff;
                el1.diff = min_diff;
                el1.it2 = el2.it1;
            }
        }
    }

   
    return 0;
}



