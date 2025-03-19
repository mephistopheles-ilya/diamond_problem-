#include <iostream>
#include <vector>
#include <algorithm>

#include <boost/program_options.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/centroid.h>
#include <CGAL/squared_distance_3.h>


#include <fenv.h>

#include "munkres.h" // from https://github.com/saebyn/munkres-cpp.git
// g++ -g -std=c++20 -fsanitize=leak,undefined,address cmp_poly.cpp -lboost_program_options -lgmp -lmunkres -L /home/ilya/munkres-cpp/build -I /home/ilya/munkres-cpp/src -o cmp_poly.out`

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;
typedef K::Point_3                                          Point_3;
typedef K::Plane_3                                          Plane;
typedef Polyhedron::Face_iterator                           Face_iterator;


Point_3 calc_centr(Face_iterator fit) {
    double all_area = 0;
    size_t num = fit->size();
    auto h = fit->halfedge();
    Point_3 center(0, 0, 0);
    Point_3 p0 = h->vertex()->point();
    h = h->next();
    Point_3 p1 = h->vertex()->point();
    h = h->next();
    Point_3 p2(0, 0, 0);
    for(size_t i = 2; i < num; ++i) {
        p2 = h->vertex()->point();
        h = h->next();
        double area = std::fabs(CGAL::squared_area(p0, p1, p2));
        Point_3 center1 = Point_3((p0.x() + p1.x() + p1.x())/3, (p0.y() + p1.y() + p2.y())/3, (p0.z() + p1.z() + p2.z())/3);
        center = Point_3(center.x() + center1.x() * area, center.y() + center1.y() * area
                , center.z() + center1.z() * area);
        all_area += area;
        p1 = p2;
    }
    return Point_3(center.x()/all_area, center.y()/all_area, center.z()/all_area);
}

Point_3 calc_norm(Face_iterator fit) {
    auto h = (*fit).halfedge();
    Plane plane  = Plane(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
    double norm_length = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
    Point_3 norm = Point_3(plane.a() / norm_length, plane.b() / norm_length, plane.c() / norm_length);
    return norm;
}

Plane calc_unit_plane(Face_iterator fit) {
    auto h = (*fit).halfedge();
    Plane plane  = Plane(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
    double norm_length = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
    Plane norm = Plane(plane.a() / norm_length, plane.b() / norm_length, plane.c() / norm_length, plane.d() / norm_length);
    return norm;
}



bool find_rundist_up_low(std::vector<Face_iterator>& planes_its
        , std::vector<Face_iterator>& rundist_planes_its
        , std::vector<Face_iterator>& up_rundist_planes_its
        , std::vector<Face_iterator>& low_rundist_planes_its
        ) {
    std::sort(planes_its.begin(), planes_its.end(), [](Face_iterator fit1, Face_iterator fit2) { 
            Plane p1 = calc_unit_plane(fit1);
            Plane p2 = calc_unit_plane(fit2);
            return p1.c() > p2.c();});

    size_t sz = planes_its.size();

    std::vector<std::pair<size_t, double>> gaps;
    gaps.reserve(sz - 1);

    for(size_t i = 0; i < sz - 1; ++i) {
        Plane p1 = calc_unit_plane(planes_its[i]);
        Plane p2 = calc_unit_plane(planes_its[i + 1]);
        double gap = p1.c() - p2.c();
        gaps.push_back({i, gap});
    }


    std::sort(gaps.begin(), gaps.end(), [](const std::pair<size_t, double>& l, const std::pair<size_t, double>& r) {
            return l.second > r.second; });
    if (gaps.size() < 2) {
        return false;
    }
    size_t begin = std::min(gaps[0].first, gaps[1].first);
    size_t before_end = std::max(gaps[0].first, gaps[1].first);
    ++begin;
    std::copy(planes_its.begin() + begin, planes_its.begin() + before_end + 1, std::back_inserter(rundist_planes_its));


    double z_start = calc_unit_plane(*rundist_planes_its.begin()).c();
    double z_end = calc_unit_plane(*std::prev(rundist_planes_its.end())).c();
    for(auto& el: planes_its) {
        Plane p = calc_unit_plane(el);
        if (p.c() > z_start) {
            up_rundist_planes_its.push_back(el);
        } else if (p.c() < z_end) {
            low_rundist_planes_its.push_back(el);
        }
    }

    return true;
}





struct face_diff {
    Face_iterator it_self;
    Face_iterator it_other;
    double diff = -1;

    double cmp_with(Face_iterator fit) const {
        Point_3 norm1 = calc_norm(it_self);
        Point_3 norm2 = calc_norm(fit);
        Point_3 center1 = calc_centr(it_self);
        Point_3 center2 = calc_centr(fit);
        return CGAL::squared_distance(norm1, norm2) + CGAL::squared_distance(center1, center2);
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


void usual_compare(std::vector<face_diff>& diffs, std::vector<Face_iterator>& poly1_its, std::vector<Face_iterator>& poly2_its) {
    for(auto& el: poly1_its) {
        face_diff fd;
        fd.it_self = el;
        diffs.push_back(fd);
    }

    for(auto& el: diffs) {
        double min_diff = std::numeric_limits<double>::max();
        double diff = 0;
        for(auto& el1: poly2_its) {
            diff = el.cmp_with(el1);
            if (diff < min_diff) {
                min_diff = diff;
                el.diff = min_diff;
                el.it_other = el1;
            }
        }
    }
}

void kuhn_compare(std::vector<face_diff>& diffs, std::vector<Face_iterator>& poly1_its, std::vector<Face_iterator>& poly2_its) {
    Matrix<double> matrix(poly1_its.size(), poly2_its.size());
    size_t nrows = matrix.rows();
    size_t ncols = matrix.columns();
    for(size_t row = 0; row < nrows; ++row) {
        face_diff fd;
        fd.it_self = poly1_its[row];
        for(size_t col = 0; col < ncols; ++col) {
            matrix(row, col) = fd.cmp_with(poly2_its[col]);
        }
    }

	Munkres<double> m;
	m.solve(matrix);

    for(auto& el: poly1_its) {
        face_diff fd;
        fd.it_self = el;
        diffs.push_back(fd);
    }

    for(size_t row = 0; row < nrows; ++row) {
        for(size_t col = 0; col < ncols; ++col) {
            if(matrix(row, col) == 0) {
                diffs[row].it_other = poly2_its[col];
                diffs[row].diff = diffs[row].cmp_with(poly2_its[col]);
                break;
            }
        }
    }
}




int main(int argc, char* argv[]) {

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

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

    std::vector<Face_iterator> poly1_all_its;
    std::vector<Face_iterator> poly1_rundits_its;
    std::vector<Face_iterator> poly1_up_rundist_its;
    std::vector<Face_iterator> poly1_low_rundist_its;

    std::vector<Face_iterator> poly2_all_its;
    std::vector<Face_iterator> poly2_rundits_its;
    std::vector<Face_iterator> poly2_up_rundist_its;
    std::vector<Face_iterator> poly2_low_rundist_its;


    for(auto fit = poly1.facets_begin(); fit != poly1.facets_end(); ++fit) {
        poly1_all_its.push_back(fit);
    }
    for(auto fit = poly2.facets_begin(); fit != poly2.facets_end(); ++fit) {
        poly2_all_its.push_back(fit);
    }

    find_rundist_up_low(poly1_all_its, poly1_rundits_its, poly1_up_rundist_its, poly1_low_rundist_its);
    find_rundist_up_low(poly2_all_its, poly2_rundits_its, poly2_up_rundist_its, poly2_low_rundist_its);

    std::vector<face_diff> diffs;
    kuhn_compare(diffs, poly1_up_rundist_its, poly2_up_rundist_its);

    std::sort(diffs.begin(), diffs.end(), [](const face_diff& fd1, const face_diff& fd2) { return fd1.diff > fd2.diff; });
    for(auto& el: diffs) {
        std::cout << el.diff << std::endl;
    }

    for(int i = 0; i < 7; ++i) {
        write_facet_ply(diffs[i].it_self, std::string("face_") + std::to_string(i) + std::string("_1.ply"));
        write_facet_ply(diffs[i].it_other, std::string("face_") + std::to_string(i) + std::string("_2.ply"));
    }

   
    return 0;
}



