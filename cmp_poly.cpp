#include <cstddef>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_set>
#include <numbers>

#include <boost/program_options.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/centroid.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <boost/property_map/vector_property_map.hpp>

#include <fenv.h>

#include "munkres.h" // from https://github.com/saebyn/munkres-cpp.git
// g++ -g -std=c++20 -fsanitize=leak,undefined,address cmp_poly.cpp -lboost_program_options -lgmp -lmunkres -L /home/ilya/munkres-cpp/build -I /home/ilya/munkres-cpp/src -o cmp_poly.out`

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;
typedef K::Point_3                                          Point_3;
typedef K::Plane_3                                          Plane;
typedef K::Vector_3                                         Vector_3;
typedef K::Segment_3                                        Segment_3;
typedef Polyhedron::Face_iterator                           Face_iterator;
typedef Polyhedron::Vertex_handle                           Vertex_handle;
typedef Polyhedron::HalfedgeDS                              HalfedgeDS;

#define _DTOR_   0.0174532925199432957692         //  PI / 180
inline double Rad2Deg(double a) { return a / _DTOR_; }
inline double Deg2Rad(double a) { return a * _DTOR_; }

//Example (1, 0, 1), (-1, 0, 1)

double Azimuth(const Point_3& v) {
    double a = std::atan2(v.x(), v.y());
    return (a >= 0) ? a : a + 2*M_PI;
}

double Slope(const Point_3 &v) {
    double d = sqrt(v.x()*v.x() + v.y()*v.y());
    return std::atan2(v.z(), d);
}


template <class HDS>
class Build_polyhedron : public CGAL::Modifier_base<HDS> {
public:
    std::vector<Face_iterator> faces;
    Build_polyhedron() {}


    void add_face(const std::vector<Face_iterator>& all_faces, size_t begin, size_t end) {
        for(size_t i = begin; i < end; ++i) {
            faces.push_back(all_faces[i]);
        }
    }
    void operator()( HDS& hds) {
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        std::unordered_map<Vertex_handle, size_t> vertices;
        size_t vertex_index = 0;
        for(auto& face: faces) {
            auto begin = face->facet_begin();
            size_t num = face->size();
            for(size_t i = 0; i < num; ++i, ++begin) {
                auto insert = vertices.insert({begin->vertex(), vertex_index}).second;
                if (insert == true) {
                    ++vertex_index;
                }
            }
        }
        B.begin_surface(vertices.size(), faces.size(), vertices.size(), 1);
        std::vector<std::pair<Vertex_handle, size_t>> v;
        std::copy(vertices.begin(), vertices.end(), std::back_inserter(v));
        std::sort(v.begin(), v.end(), [](const std::pair<Vertex_handle, size_t>& a
                    , const std::pair<Vertex_handle, size_t>& b) { return a.second < b.second;});
        for(const auto& el: v) {
            B.add_vertex(el.first->point());
        }
        for(auto& face: faces) {
            B.begin_facet();
            auto begin = face->facet_begin();
            size_t num = face->size();
            for(size_t i = 0; i < num; ++i, ++begin) {
                B.add_vertex_to_facet(vertices.find(begin->vertex())->second);
            }
            B.end_facet();
        }
        B.end_surface();
    }
};

void print_faces_as_polyhedron_ply(std::vector<Face_iterator>& all_faces, size_t begin, size_t end, std::string filename
        , double precison = 12.) {
    //std::cout << "HERE " << all_faces.size()  << ' ' << filename << std::endl;
    Polyhedron polyhedron;
    Build_polyhedron<HalfedgeDS> build_poly;
    build_poly.add_face(all_faces, begin, end);
    polyhedron.delegate(build_poly);

    std::ofstream of(filename);
    of << std::fixed << std::setprecision(precison);
    CGAL::IO::write_PLY(of, polyhedron);
}


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
        Point_3 center1 = Point_3((p0.x() + p1.x() + p2.x())/3, (p0.y() + p1.y() + p2.y())/3, (p0.z() + p1.z() + p2.z())/3);
        center = Point_3(center.x() + center1.x() * area, center.y() + center1.y() * area
                , center.z() + center1.z() * area);
        all_area += area;
        p1 = p2;
    }
    return Point_3(center.x()/all_area, center.y()/all_area, center.z()/all_area);
}


#if 0
Point_3 calc_norm(Face_iterator fit) {
    static std::vector<Point_3> points(100);

    if (fit->size() > points.capacity()) {
        points.resize(fit->size());
    }
    points.clear();

    auto h = (*fit).halfedge();
    size_t pn = fit->size();
    size_t i = 0;
    for(; i < pn; ++i, h = h->next()) {
        points.push_back(h->vertex()->point());
    }
    Plane best_plane;
    Point_3 centroid;

    CGAL::linear_least_squares_fitting_3( points.begin(), points.end(), best_plane, centroid, CGAL::Dimension_tag<0>());
    Vector_3 unit_normal = best_plane.orthogonal_vector();
    Vector_3 oriented_normal = Plane(points[0], points[1], points[3]).orthogonal_vector();
    double prod = CGAL::scalar_product(unit_normal, oriented_normal);
    if (prod < 0) {
        unit_normal = unit_normal * (-1);
    }
    unit_normal = unit_normal / unit_normal.squared_length();
    return Point_3(unit_normal.x(), unit_normal.y(), unit_normal.z());
}
#endif

#if 0
Point_3 calc_norm(Face_iterator fit) {
    static std::vector<Segment_3> segments(100);

    if (fit->size() > segments.capacity()) {
        segments.resize(fit->size());
    }
    segments.clear();

    auto h = (*fit).halfedge();
    size_t pn = fit->size();
    size_t i = 0;
    for(; i < pn; ++i, h = h->next()) {
        segments.push_back(Segment_3(h->vertex()->point(), h->next()->vertex()->point()));
    }
    Plane plane;
    CGAL::linear_least_squares_fitting_3(segments.begin(), segments.end(), plane, CGAL::Dimension_tag<1>());
    Vector_3 v = Vector_3(plane.a(), plane.b(), plane.c());
    double length = std::sqrt(v.squared_length());
    v =  v / length;
    Plane or_plane = Plane(segments[0][0], segments[0][1], segments[1][1]);
    Vector_3  or_v = Vector_3(or_plane.a(), or_plane.b(), or_plane.c());
    length = std::sqrt(or_v.squared_length());
    or_v = or_v /  length;
    double prod = CGAL::scalar_product(v, or_v);
    if (prod < 0) {
        v =  -v;
    }

    return Point_3(v.x(), v.y(), v.z());

}
#endif


#if 1
Point_3 calc_norm(Face_iterator fit) {
    static std::vector<Point_3> points(100);

    if (fit->size() > points.capacity()) {
        points.resize(fit->size());
    }
    points.clear();

    auto h = (*fit).halfedge();
    size_t pn = fit->size();
    size_t i = 0;
    for(; i < pn; ++i, h = h->next()) {
        points.push_back(h->vertex()->point());
    }
    //Plane plane;
    //CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), plane, CGAL::Dimension_tag<0>());

    Vector_3 unit_normal(0, 0, 0);

    Point_3 a(0, 0, 0);
    Point_3 b(0, 0, 0);
    Point_3 c(0, 0, 0);

    double max_dist = -1;
    double current_dist = -1;
    size_t index_a = 0;
    size_t index_b = 0;
    for(i = 0; i < pn - 1; ++i) {
        current_dist = CGAL::squared_distance(points[i], points[i + 1]);
        if (current_dist > max_dist) {
            a = points[i];
            b = points[i + 1];
            index_a = i;
            index_b = i + 1;
            max_dist = current_dist;
        }
    }
    current_dist = CGAL::squared_distance(points[pn - 1], points[0]);
    if (current_dist > max_dist) {
        a = points[pn - 1];
        b = points[0];
        index_a = pn - 1;
        index_b = 0;
    }

    max_dist = -1;
    K::Line_3 line(a, b);

    for(i = 0; i < pn; ++i) {
        if (i == index_a || i == index_b) {
            continue;
        }
        current_dist = CGAL::squared_distance(points[i], line);
        if (current_dist > max_dist) {
            c = points[i];
            max_dist = current_dist;
        }
    }

    unit_normal = CGAL::unit_normal(a, b, c);

    //Point_3 centr = CGAL::centroid(a, b, c);
    //Vector_3 v_centr = Vector_3(Point_3(0, 0, 0), centr);
    //double prod = CGAL::scalar_product(v_centr, unit_normal);
#if 0
    Plane plane = Plane(points[0], points[1], points[2]);
    double norm_length = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
    Plane norm_plane = Plane(plane.a() / norm_length, plane.b() / norm_length, plane.c() / norm_length, plane.d() / norm_length);
    Vector_3 unit_v = Vector_3(norm_plane.a(), norm_plane.b(), norm_plane.c());
    double prod = CGAL::scalar_product(unit_v, unit_normal);
    if (prod < 0) {
        unit_normal = -unit_normal;
    }
#endif

    return Point_3(unit_normal.x(), unit_normal.y(), unit_normal.z());
}
#endif


#if 0
Plane calc_unit_plane(Face_iterator fit) {
    static std::vector<Point_3> points(100);

    if (fit->size() > points.capacity()) {
        points.resize(fit->size());
    }
    points.clear();

    auto h = (*fit).halfedge();
    size_t pn = fit->size();
    size_t i = 0;
    for(; i < pn; ++i, h = h->next()) {
        points.push_back(h->vertex()->point());
    }


    Point_3 a(0, 0, 0);
    Point_3 b(0, 0, 0);
    Point_3 c(0, 0, 0);

    double max_dist = -1;
    double current_dist = -1;
    size_t index_a = 0;
    size_t index_b = 0;
    for(i = 0; i < pn - 1; ++i) {
        current_dist = CGAL::squared_distance(points[i], points[i + 1]);
        if (current_dist > max_dist) {
            a = points[i];
            b = points[i + 1];
            index_a = i;
            index_b = i + 1;
            max_dist = current_dist;
        }
    }
    current_dist = CGAL::squared_distance(points[pn - 1], points[0]);
    if (current_dist > max_dist) {
        a = points[pn - 1];
        b = points[0];
        index_a = pn - 1;
        index_b = 0;
    }

    max_dist = -1;
    K::Line_3 line(a, b);

    for(i = 0; i < pn; ++i) {
        if (i == index_a || i == index_b) {
            continue;
        }
        current_dist = CGAL::squared_distance(points[i], line);
        if (current_dist > max_dist) {
            c = points[i];
            max_dist = current_dist;
        }
    }

    Plane plane = Plane(a, b, c);
    double norm_length = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
    Plane norm = Plane(plane.a() / norm_length, plane.b() / norm_length, plane.c() / norm_length, plane.d() / norm_length);
    return norm;
}
#endif



bool find_rundist_up_low(std::vector<Face_iterator>& planes_its
        , std::vector<Face_iterator>& rundist_planes_its
        , std::vector<Face_iterator>& up_rundist_planes_its
        , std::vector<Face_iterator>& low_rundist_planes_its
        ) {
    std::sort(planes_its.begin(), planes_its.end(), [](Face_iterator fit1, Face_iterator fit2) {
            Point_3 p1 = calc_norm(fit1);
            Point_3 p2 = calc_norm(fit2);
            return p1.z() > p2.z();});
#if 0
    for(auto& pit : planes_its) {
        fprintf(stderr,"%.10lf\n",  calc_norm(pit).z());
    }
    std::cerr << std::endl;
#endif

    size_t sz = planes_its.size();

    std::vector<std::pair<size_t, double>> gaps;
    gaps.reserve(sz - 1);

    for(size_t i = 0; i < sz - 1; ++i) {
        Point_3 p1 = calc_norm(planes_its[i]);
        Point_3 p2 = calc_norm(planes_its[i + 1]);
        double gap = p1.z() - p2.z();
        gaps.push_back({i, gap});
    }

    size_t end_of_the_top = 0;
    for(size_t i = 0; i < sz - 2; ++i) {
        if (gaps[i + 1].second < gaps[i].second) {
            end_of_the_top = i;
            break;
        }
    }
    if (calc_norm(planes_its[end_of_the_top]).z() < 0.9) {
        end_of_the_top = -1;
    }
#if 0
    static int counter = 0;
    std::cerr << end_of_the_top << std::endl;
    std::vector<Face_iterator> top;
    for(size_t i = 0; i <= end_of_the_top; ++i) {
        top.push_back(planes_its[i]);
    }
    print_faces_as_polyhedron_ply(top, 0, top.size(), std::string("top_") + std::to_string(counter) + std::string(".ply"));
    ++counter;
#endif


    std::sort(gaps.begin(), gaps.end(), [](const std::pair<size_t, double>& l, const std::pair<size_t, double>& r) {
            return l.second > r.second; });
    if (gaps.size() < 2) {
        return false;
    }
#if 1
    size_t max_gap1 = 0,  max_gap2 = 0;
    for(size_t i = 0; i < gaps.size(); ++i) {
        if (std::fabs(calc_norm(planes_its[gaps[i].first]).z()) < 0.7 || std::fabs(calc_norm(planes_its[gaps[i].first + 1]).z()) < 0.7) {
            max_gap1 = i;
            break;
        }
    }
    for(size_t i = max_gap1 + 1; i < gaps.size(); ++i) {
        if (std::fabs(calc_norm(planes_its[gaps[i].first]).z()) < 0.7 || std::fabs(calc_norm(planes_its[gaps[i].first + 1]).z()) < 0.7) {
            max_gap2 = i;
            break;
        }
    }
#endif

    size_t begin = std::min(gaps[max_gap1].first, gaps[max_gap2].first);
    size_t before_end = std::max(gaps[max_gap1].first, gaps[max_gap2].first);
    ++begin;
    std::copy(planes_its.begin() + begin, planes_its.begin() + before_end + 1, std::back_inserter(rundist_planes_its));

    //std::cout << calc_norm(planes_its[begin]).z() << ' ' << calc_norm(planes_its[before_end]).z() << std::endl;

    double z_start = calc_norm(*rundist_planes_its.begin()).z();
    double z_end = calc_norm(*std::prev(rundist_planes_its.end())).z();
    for(size_t i = end_of_the_top + 1; i < planes_its.size(); ++i) {
        auto& el = planes_its[i];
        Point_3 p = calc_norm(el);
        if (p.z() > z_start) {
            up_rundist_planes_its.push_back(el);
        } else if (p.z() < z_end) {
            low_rundist_planes_its.push_back(el);
        }
    }
    //std::cout << "UP sz = " << up_rundist_planes_its.size() << ' ';
    //std::cout << "R sz = " << rundist_planes_its.size() << ' ';
    //std::cout << "LOW sz = " << low_rundist_planes_its.size() << std::endl;

    return true;
}





struct face_diff {
    Face_iterator it_self;
    Face_iterator it_other;
    double diff_2 = -1;
    double diff_slope = -1;
    double diff_azimuth = -1;
    double diff_angle = -1;
    bool matched = false;

    double cmp_with(Face_iterator fit) const {
        Point_3 norm1 = calc_norm(it_self);
        Point_3 norm2 = calc_norm(fit);
        Point_3 center1 = calc_centr(it_self);
        Point_3 center2 = calc_centr(fit);
        return CGAL::squared_distance(norm1, norm2) + CGAL::squared_distance(center1, center2);
    }
    double cmp_normal_2(Face_iterator fit) const {
        Point_3 norm1 = calc_norm(it_self);
        Point_3 norm2 = calc_norm(fit);
        return CGAL::squared_distance(norm1, norm2);
    }
    double cmp_normal_azimuth(Face_iterator fit) const {
        const double pi = std::numbers::pi;
        Point_3 norm1 = calc_norm(it_self);
        Point_3 norm2 = calc_norm(fit);
        //std::cout << norm1.x() << ' ' << norm1.y() << ' ' << norm1.z() << std::endl;
        //std::cout << norm2.x() << ' ' << norm2.y() << ' ' << norm2.z() << std::endl;
        if (std::fabs(norm1.x()) < 1e-2 || std::fabs(norm1.y()) < 1e-2 || std::fabs(norm2.x()) < 1e-2 || std::fabs(norm2.y()) < 1e-2) {
            return 0;
        }
        double azimuth1 = Azimuth(norm1);
        double azimuth2 = Azimuth(norm2);
        //printf("%lf %lf\n", Rad2Deg(azimuth1), Rad2Deg(azimuth2));
        double d = std::fabs(azimuth1 - azimuth2);
        if (2 * pi - d < d) {
            d = 2 * pi - d;
        }
        //if (d < 0) {
        //    std::cout << "Les than zero" << std::endl;
        //}
        return Rad2Deg(d);
    }
    double cmp_normal_slope(Face_iterator fit) const {
        Point_3 norm1 = calc_norm(it_self);
        Point_3 norm2 = calc_norm(fit);
        double slope1 = Slope(norm1);
        double slope2 = Slope(norm2);
        double d = std::fabs(slope1 - slope2);
        return Rad2Deg(d);
    }
    double cmp_angle(Face_iterator fit) const {
        Point_3 norm1 = calc_norm(it_self);
        Point_3 norm2 = calc_norm(fit);
        Vector_3 v1 = Vector_3(norm1.x(), norm1.y(), norm1.z());
        Vector_3 v2 = Vector_3(norm2.x(), norm2.y(), norm2.z());
        double v1_norm = CGAL::sqrt(v1.squared_length());
        double v2_norm = CGAL::sqrt(v2.squared_length());
        double cos_theta = (v1 * v2) / (v1_norm * v2_norm);
        cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
        double angle = std::acos(cos_theta);
        return Rad2Deg(angle);
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
                el.diff_2 = min_diff;
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
            if(std::fabs(matrix(row, col)) < 1e-16) {
                diffs[row].it_other = poly2_its[col];
                diffs[row].diff_2 = diffs[row].cmp_normal_2(poly2_its[col]);
                diffs[row].diff_slope = diffs[row].cmp_normal_slope(poly2_its[col]);
                diffs[row].diff_azimuth = diffs[row].cmp_normal_azimuth(poly2_its[col]);
                diffs[row].diff_angle = diffs[row].cmp_angle(poly2_its[col]);
                diffs[row].matched = true;
                break;
            }
        }
    }
#if 0
    for(size_t row = 0; row < nrows; ++row) {
        for(size_t col = 0; col < ncols; ++col) {
            fprintf(stderr, "%e ", matrix(row, col));
        }
        fprintf(stderr, "\n");
    }
#endif
}


bool check_match(const std::vector<face_diff>& diffs) {
    for(const auto& el: diffs) {
        if (el.matched == false) {
            return false;
        }
    }
    return true;
}

double mean_l2_value(const std::vector<face_diff>& diffs, double face_diff::*diff) {
    double val = 0;
    int count = 0;
    for(const auto& el: diffs) {
        double diff_value = std::invoke(diff, el);
        if (diff_value >= 0) {
            ++count;
            val +=  diff_value * diff_value;
        }
    }
    val /= count;
    val = std::sqrt(val);
    return val;
}
double mean_value(const std::vector<face_diff>& diffs, double face_diff::*diff) {
    double val = 0;
    int count = 0;
    for(const auto& el: diffs) {
        double diff_value = std::invoke(diff, el);
        if (diff_value >= 0) {
            ++count;
            val += diff_value;
        }
    }
    assert(count > 0 && "Number of elements must be greater tha zero");
    val /= count;
    return val;
}

double std_deviation(const std::vector<face_diff>& diffs, double mean_val, double face_diff::*diff) {
    double val = 0;
    int count = 0;
    for(const auto& el: diffs) {
        double diff_value = std::invoke(diff, el);
        if (diff_value >= 0) {
            ++count;
            val += (diff_value - mean_val) * (diff_value - mean_val);
        }
    }
    assert(count > 0 && "Number of elements must be greater tha zero");
    val /= count;
    return std::sqrt(val);
}


int get_face_number(Polyhedron& poly, Face_iterator it) {
    int counter = 0;
    for(auto fit = poly.facets_begin(); fit != poly.facets_end(); ++fit, ++counter) {
        if (it == fit) {
            return counter;
        }
    }
    return -1;
}

void print_results(Polyhedron& poly1, Polyhedron& poly2, std::vector<face_diff>& diffs, std::string part, double face_diff::*diff
        , std::string file1, std::string file2, int edges_d, int face_d) {
    std::sort(diffs.begin(), diffs.end(), [diff](const face_diff& fd1, const face_diff& fd2) { return std::invoke(diff, fd1) 
            > std::invoke(diff, fd2); });
    std::string name1 = file1.substr(file1.find("Reflect"));
    std::string name2 = file2.substr(file2.find("Reflect"));
    std::cout << name1 << " /  " << name2 << std::endl;
    if (check_match(diffs) == false) {
        std::cout << "Not all faces from " << part << " are matched" << std::endl;
    } else {
        std::cout << "All faces from " << part << " are matched" << std::endl;
    }
    std::cout << "Edges : " << edges_d << std::endl;
    std::cout << "Faces : " << face_d << std::endl;
    if (diff == &face_diff::diff_azimuth) {
        std::cout << "Azimuth in degrees" << std::endl;
    }
    if (diff == &face_diff::diff_slope) {
        std::cout << "Slope in degrees" << std::endl;
    }
    if (diff == &face_diff::diff_2) {
        std::cout << "Distance in microns" << std::endl;
    }
    if (diff == &face_diff::diff_angle) {
        std::cout << "Angle in degrees" << std::endl;
    }

    double mean = mean_value(diffs, diff);
    double max = std::invoke(diff, diffs[0]);
    double std_dev = std_deviation(diffs, mean, diff);
    std::cout << part << " max = " << max << std::endl;
    std::cout << part << " mean = " << mean << std::endl;
    std::cout << part << " std.dev. = " << std_dev << std::endl;

    std::cout << "Peaks in " << part << std::endl;
    for(size_t i = 0; i < std::min(diffs.size(), 7ul); ++i) {
        int num1 = get_face_number(poly1, diffs[i].it_self);
        int num2 = get_face_number(poly2, diffs[i].it_other);
        std::cout << "Diff= " << std::invoke(diff, diffs[i]) << " ";
        std::cout << "for face : ";
        std::cout << num1 << " / " << num2 << std::endl;
    }
    std::cout << std::endl;

#if 0
    for(auto& el: diffs) {
        std::cout << std::invoke(diff, el) << std::endl;
    }
    std::cout << std::endl;
#endif
}

int calculate_num_of_edges(std::vector<Face_iterator>& poly_faces, std::unordered_set<Polyhedron::Vertex_handle>& rundist_vertices) {
    int edges_in_poly = 0;
    for(auto& face: poly_faces) {
        auto h = face->halfedge();
        size_t nh = face->size();
        int num_hedges = 0;
        for(size_t i = 0; i < nh; ++i, h = h->next()) {
            auto v1 = h->vertex();
            auto v2 = h->opposite()->vertex();
            auto end = rundist_vertices.end();
            if (rundist_vertices.find(v1) != end && rundist_vertices.find(v2) != end) {
                continue;
            }
            num_hedges++;
        }
        edges_in_poly += num_hedges;
    }
    edges_in_poly /= 2;
    return edges_in_poly;
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
    std::vector<Face_iterator> poly1_rundist_its;
    std::vector<Face_iterator> poly1_up_rundist_its;
    std::vector<Face_iterator> poly1_low_rundist_its;

    std::vector<Face_iterator> poly2_all_its;
    std::vector<Face_iterator> poly2_rundist_its;
    std::vector<Face_iterator> poly2_up_rundist_its;
    std::vector<Face_iterator> poly2_low_rundist_its;


    for(auto fit = poly1.facets_begin(); fit != poly1.facets_end(); ++fit) {
        poly1_all_its.push_back(fit);
    }
    for(auto fit = poly2.facets_begin(); fit != poly2.facets_end(); ++fit) {
        poly2_all_its.push_back(fit);
    }

    find_rundist_up_low(poly1_all_its, poly1_rundist_its, poly1_up_rundist_its, poly1_low_rundist_its);
    std::unordered_set<Polyhedron::Vertex_handle> rundist_vertices_poly1;
    for(auto& el: poly1_rundist_its) {
        auto begin = el->facet_begin();
        size_t sz = el->size();
        for(size_t i = 0; i < sz; ++i, ++begin) {
            rundist_vertices_poly1.insert(begin->vertex());
        }
    }
    int edges_in_poly1_up_rundist = 0;
    int edges_in_poly1_low_rundist = 0;
    edges_in_poly1_up_rundist = calculate_num_of_edges(poly1_up_rundist_its, rundist_vertices_poly1);
    edges_in_poly1_low_rundist = calculate_num_of_edges(poly1_low_rundist_its, rundist_vertices_poly1);

    find_rundist_up_low(poly2_all_its, poly2_rundist_its, poly2_up_rundist_its, poly2_low_rundist_its);
    std::unordered_set<Polyhedron::Vertex_handle> rundist_vertices_poly2;
    for(auto& el: poly2_rundist_its) {
        auto begin = el->facet_begin();
        size_t sz = el->size();
        for(size_t i = 0; i < sz; ++i, ++begin) {
            rundist_vertices_poly2.insert(begin->vertex());
        }
    }
    int edges_in_poly2_up_rundist = 0;
    int edges_in_poly2_low_rundist = 0;
    edges_in_poly2_up_rundist = calculate_num_of_edges(poly2_up_rundist_its, rundist_vertices_poly2);
    edges_in_poly2_low_rundist = calculate_num_of_edges(poly2_low_rundist_its, rundist_vertices_poly2);

#if 1
    print_faces_as_polyhedron_ply(poly1_rundist_its, 0, poly1_rundist_its.size(), file1.substr(file1.find("Reflect")) + std::string("_rundist.ply"));
    print_faces_as_polyhedron_ply(poly1_up_rundist_its, 0, poly1_up_rundist_its.size(), file1.substr(file1.find("Reflect")) + std::string("_up_rundist.ply"));
    print_faces_as_polyhedron_ply(poly1_low_rundist_its, 0, poly1_low_rundist_its.size(), file1.substr(file1.find("Reflect")) + std::string("_low_rundist.ply"));

    print_faces_as_polyhedron_ply(poly2_rundist_its, 0, poly2_rundist_its.size(), file2.substr(file2.find("Reflect")) + std::string("_rundist.ply"));
    print_faces_as_polyhedron_ply(poly2_up_rundist_its, 0, poly2_up_rundist_its.size(), file2.substr(file2.find("Reflect")) + std::string("_up_rundist.ply"));
    print_faces_as_polyhedron_ply(poly2_low_rundist_its, 0, poly2_low_rundist_its.size(), file2.substr(file2.find("Reflect")) + std::string("_low_rundist.ply"));
#endif


    std::vector<face_diff> diffs_up_rundist;
    kuhn_compare(diffs_up_rundist, poly1_up_rundist_its, poly2_up_rundist_its);

    std::vector<face_diff> diffs_low_rundist;
    kuhn_compare(diffs_low_rundist, poly1_low_rundist_its, poly2_low_rundist_its);

#if 0
    print_results(poly1, poly2, diffs_up_rundist, "pavilion", &face_diff::diff_azimuth, file1, file2
            , edges_in_poly1_up_rundist - edges_in_poly2_up_rundist
            , static_cast<int>(poly1_up_rundist_its.size()) - static_cast<int>(poly2_up_rundist_its.size()));
#endif
#if 0
    print_results(poly1, poly2, diffs_up_rundist, "pavilion", &face_diff::diff_slope, file1, file2
            , edges_in_poly1_up_rundist - edges_in_poly2_up_rundist
            , static_cast<int>(poly1_up_rundist_its.size()) - static_cast<int>(poly2_up_rundist_its.size()));
#endif
#if 0
    print_results(poly1, poly2, diffs_up_rundist, "pavilion", &face_diff::diff_angle, file1, file2
            , edges_in_poly1_up_rundist - edges_in_poly2_up_rundist
            , static_cast<int>(poly1_up_rundist_its.size()) - static_cast<int>(poly2_up_rundist_its.size()));
#endif
#if 0
    print_results(poly1, poly2, diffs_low_rundist, "crown", &face_diff::diff_azimuth, file1, file2
            , edges_in_poly1_low_rundist - edges_in_poly2_low_rundist
            , static_cast<int>(poly1_low_rundist_its.size()) - static_cast<int>(poly2_low_rundist_its.size()));
#endif
#if 0
    print_results(poly1, poly2, diffs_low_rundist, "crown", &face_diff::diff_slope, file1, file2
            , edges_in_poly1_low_rundist - edges_in_poly2_low_rundist
            , static_cast<int>(poly1_low_rundist_its.size()) - static_cast<int>(poly2_low_rundist_its.size()));
#endif
#if 1
    print_results(poly1, poly2, diffs_low_rundist, "crown", &face_diff::diff_angle, file1, file2
            , edges_in_poly1_low_rundist - edges_in_poly2_low_rundist
            , static_cast<int>(poly1_low_rundist_its.size()) - static_cast<int>(poly2_low_rundist_its.size()));
#endif



#define part diffs_low_rundist
    for(int i = 0; i < std::min(7ul, part.size()); ++i) {
        write_facet_ply(part[i].it_self, std::string("face_") + std::to_string(i) + std::string("_1.ply"));
        write_facet_ply(part[i].it_other, std::string("face_") + std::to_string(i) + std::string("_2.ply"));
#if 1
        Point_3 norm1 = calc_norm(part[i].it_self);
        Point_3 norm2 = calc_norm(part[i].it_other);
        face_diff fd;
        fd.it_self = part[i].it_self;
        std::cout << "????   " << fd.cmp_angle(part[i].it_other) << std::endl;
        std::cout << get_face_number(poly1, part[i].it_self) << " / " << get_face_number(poly2, part[i].it_other) << std::endl; 
        Point_3 center1 = calc_centr(part[i].it_self);
        Point_3 center2 = calc_centr(part[i].it_other);

        Point_3 p(center1.x() + norm1.x(), center1.y() + norm1.y(), center1.z() + norm1.z());
        write_segment_ply(center1, p ,std::string("norm_") + std::to_string(i) + std::string("_1.ply"));

        Point_3 pp(center2.x() + norm2.x(), center2.y() + norm2.y(), center2.z() + norm2.z());

        write_segment_ply(center2, pp, std::string("norm_") + std::to_string(i) + std::string("_2.ply"));
#endif
    }



    return 0;
}



