#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>

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
    Point_3 norm1;
    Point_3 center1;
    typename Polyhedron::Face_iterator it1;
    typename Polyhedron::Face_iterator it2;
    Point_3 norm2;
    Point_3 center2;
    double diff = -1;

    double operator()(const norm_center& nc) const {
        return CGAL::squared_distance(nc.norm1, norm1) + CGAL::squared_distance(nc.center1, center1);
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

void write_point_ply(Point_3 p, const std::string& filename) {
    std::ofstream out(filename);

    out << "ply" << std::endl;
    out << "format ascii 1.0" << std::endl;
    out << "element vertex 1" << std::endl;
    out << "property double x" << std::endl;
    out << "property double y" << std::endl;
    out << "property double z" << std::endl;
    out << "end_header" << std::endl;

    out << p.x() << " " << p.y() << " " << p.z() << std::endl;

    out.close();
}

void write_segment_ply(Point_3 p1, Point_3 p2, const std::string& filename) {
    std::ofstream out(filename);

    out << "ply" << std::endl;
    out << "format ascii 1.0" << std::endl;
    out << "element vertex " << 2 << std::endl;
    out << "property double x" << std::endl;
    out << "property double y" << std::endl;
    out << "property double z" << std::endl;
    out << "element edge " << 1 << std::endl;
    out << "property int vertex1" << std::endl;
    out << "property int vertex2" << std::endl;
    out << "end_header" << std::endl;

    out << p1.x() << " " << p1.y() << " " << p1.z() << std::endl;
    out << p2.x() << " " << p2.y() << " " << p2.z() << std::endl;

    out << 0 << " " << 1 << std::endl;

    out.close();
}

Point_3 calc_centr(std::vector<Point_3>& points) {
    size_t num = points.size();
    Point_3 center(0, 0, 0);
    double all_area = 0;
    for(size_t i = 2; i < num; ++i) {
        double area = std::fabs(CGAL::squared_area(points[0], points[i-1], points[i]));
        Point_3 center1 = Point_3((points[0].x() + points[i-1].x() + points[i].x())/3
                , (points[0].y() + points[i-1].y() + points[i].y())/3
                , (points[0].z() + points[i-1].z() + points[i].z())/3);
        center = Point_3(center.x() + center1.x() * area, center.y() + center1.y() * area
                , center.z() + center1.z() * area);
        all_area += area;
    }
    return Point_3(center.x()/all_area, center.y()/all_area, center.z()/all_area);
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
        std::vector<Point_3> p;
        auto begin = boost::make_transform_iterator(fit->facet_begin(), to_point());
        std::copy_n(begin, fit->size(), std::back_inserter(p));
        //auto end_ = boost::make_transform_iterator(fit->facet_begin(), to_point());
        //std::advance(end, fit->size() - 1);
#if 0
        Point_3 center(0, 0, 0);
        size_t amount = fit->size();
        for(size_t i = 0; i < amount; ++i, ++begin) {
            center = Point_3(center.x() + (*begin).x(), center.y() + (*begin).y(), center.z() + (*begin).z());
        }
        center = Point_3(center.x()/amount, center.y()/amount, center.z()/amount);
#endif
    
        //Point_3 center = CGAL::centroid(p.begin(), p.end());
        Point_3 center = calc_centr(p);
        auto h = (*fit).halfedge();
        Plane plane  = Plane(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
        double norm_length = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
        Point_3 norm = Point_3(plane.a() / norm_length, plane.b() / norm_length, plane.c() / norm_length);
        mer1.emplace_back(norm, center, fit);
    }
    std::vector<norm_center> mer2;
    for(auto fit = poly2.facets_begin(); fit != poly2.facets_end(); ++fit) {
        std::vector<Point_3> p;
        auto begin = boost::make_transform_iterator(fit->facet_begin(), to_point());
        std::copy_n(begin, fit->size(), std::back_inserter(p));;
        //auto end = boost::make_transform_iterator(fit->facet_begin(), to_point());
        //std::advance(end, fit->size() - 1);
#if 0
        Point_3 center(0, 0, 0);
        size_t amount = fit->size();
        for(size_t i = 0; i < amount; ++i, ++begin) {
            center = Point_3(center.x() + (*begin).x(), center.y() + (*begin).y(), center.z() + (*begin).z());
        }
        center = Point_3(center.x()/amount, center.y()/amount, center.z()/amount);
#endif
   
        //Point_3 center = CGAL::centroid(p.begin(), p.end());
        Point_3 center = calc_centr(p);
        auto h = (*fit).halfedge();
        Plane plane  = Plane(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
        double norm_length = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
        Point_3 norm = Point_3(plane.a() / norm_length, plane.b() / norm_length, plane.c() / norm_length);
        mer2.emplace_back(norm, center, fit);
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
                el1.center2 = el2.center1;
                el1.norm2 = el2.norm1;
            }
        }
    }
    std::sort(mer1.begin(), mer1.end(), [](const norm_center& c1, const norm_center& c2) { return c1.diff > c2.diff; });
    for(auto& el: mer1) {
        std::cout << el.diff << std::endl;
    }

    for(int i = 0; i < 1; ++i) {
        write_facet_ply(mer1[i].it1, std::string("face_") + std::to_string(i) + std::string("_1.ply"));
        write_facet_ply(mer1[i].it2, std::string("face_") + std::to_string(i) + std::string("_2.ply"));
        //write_point_ply(mer1[i].center1, std::string("centr_") + std::to_string(i) + std::string("_1.ply"));
        //write_point_ply(mer1[i].center2, std::string("centr_") + std::to_string(i) + std::string("_2.ply"));
        Point_3 p(mer1[i].center1.x() + mer1[i].norm1.x(), mer1[i].center1.y() + mer1[i].norm1.y()
                , mer1[i].center1.z() + mer1[i].norm1.z());
        write_segment_ply(mer1[i].center1, p
                ,std::string("norm_") + std::to_string(i) + std::string("_1.ply"));
        Point_3 pp(mer1[i].center2.x() + mer1[i].norm2.x(), mer1[i].center2.y() + mer1[i].norm2.y()
                , mer1[i].center2.z() + mer1[i].norm2.z());
        write_segment_ply(mer1[i].center2, pp
                ,std::string("norm_") + std::to_string(i) + std::string("_2.ply"));

    }

   
    return 0;
}



