#include <iostream>
#include <ostream>
#include <vector>
#include <map>
#include <limits>
#include <tuple>
#include <unordered_set>
#include <random>

#include <boost/program_options.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <fenv.h>

#define DEBUG


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Exact_predicates_exact_constructions_kernel   K;
//typedef CGAL::Simple_cartesian<double>                      K;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;
typedef K::Point_3                                          Point_3;
typedef K::Plane_3                                          Plane;
typedef K::Vector_3                                         Vector_3;
typedef Polyhedron::Face_iterator                           Face_iterator; 
typedef Polyhedron::Face_handle                             Face_handle;
typedef Polyhedron::Vertex_handle                           Vertex_handle;
typedef Polyhedron::HalfedgeDS                              HalfedgeDS;


template <class HDS>
class Build_polyhedron : public CGAL::Modifier_base<HDS> {
public:
    std::vector<Face_iterator> faces;
    Build_polyhedron() {}

    template <typename Smt>
    void add_face(const std::vector<std::pair<Smt, Face_iterator>>& all_faces, size_t begin, size_t end) {
        for(size_t i = begin; i < end; ++i) {
            faces.push_back(all_faces[i].second);
        }
    }
    void operator()( HDS& hds) {
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        std::map<Vertex_handle, size_t> vertex_map;
        size_t vertex_index = 0;
        for(auto& face: faces) {
            auto begin = face->facet_begin();
            size_t num = face->size();
            for(size_t i = 0; i < num; ++i, ++begin) {
                if (vertex_map.find(begin->vertex()) == vertex_map.end()) {
                    vertex_map[begin->vertex()] = vertex_index;
                    ++vertex_index;
                }
            }
        }
        B.begin_surface(vertex_map.size(), faces.size(), vertex_map.size(), 1);
        std::vector<std::pair<Vertex_handle, size_t>> v;
        std::copy(vertex_map.begin(), vertex_map.end(), std::back_inserter(v));
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
                B.add_vertex_to_facet(vertex_map[begin->vertex()]);
            }
            B.end_facet();
        }
        B.end_surface();
    }
};

template <typename Smt>
void print_faces_as_polyhedron_ply(std::vector<std::pair<Smt, Face_iterator>>& all_faces, size_t begin, size_t end, std::string filename
        , double precison = 12.) {
    Polyhedron polyhedron;
    Build_polyhedron<HalfedgeDS> build_poly;
    build_poly.add_face(all_faces, begin, end);
    polyhedron.delegate(build_poly);

    std::ofstream of(filename);
    of << std::fixed << std::setprecision(precison);
    CGAL::IO::write_PLY(of, polyhedron);
}
void print_polyhedron_ply(Polyhedron& polyhedron, std::string filename, double precision=12.) {
        std::ofstream of(filename);
        of << std::fixed << std::setprecision(precision);
        CGAL::IO::write_PLY(of, polyhedron);
}




bool read_planes_from_file(std::vector<Plane>& planes, std::string filename
        , std::string start = "indices_of_vertices"
        , std::string stop = "#") {
    std::ifstream in(filename);
    if (in.is_open() == false) {
        std::cerr << "Canot open file: " << filename << std::endl;
        return false;
    }
    std::string line;
    for(;;) {
        std::getline(in, line);
        if (line.find(start) != std::string::npos) 
            break;
    }
    std::getline(in, line);

    while(line.find(stop) == std::string::npos) {
        double a, b, c, d;
        int unused1, unused2;
        sscanf(line.c_str(), "%d %d %lf %lf %lf %lf", &unused1, &unused2, &a, &b, &c, &d);
        Plane plane(a, b, c, d);
        planes.push_back(plane);
        std::getline(in, line);
        std::getline(in, line);
    }
    in.close();
    return true;
}

void find_rundist_up_low(std::vector<std::pair<Plane, Face_iterator>>& planes_its
        , std::vector<std::pair<Plane, Face_iterator>>& rundist_planes_its
        , std::vector<std::pair<Plane, Face_iterator>>& up_rundist_planes_its
        , std::vector<std::pair<Plane, Face_iterator>>& low_rundist_planes_its
        ) {
    std::sort(planes_its.begin(), planes_its.end(), [](const std::pair<Plane, Face_iterator>& p1
                , const std::pair<Plane, Face_iterator>& p2)
            { return p1.first.c() > p2.first.c();});

    size_t sz = planes_its.size();

    std::vector<std::pair<size_t, double>> gaps;
    gaps.reserve(sz - 1);

    for(size_t i = 0; i < sz - 1; ++i) {
        double gap = planes_its[i].first.c() - planes_its[i + 1].first.c();
        gaps.push_back({i, gap});
    }


    std::sort(gaps.begin(), gaps.end(), [](const std::pair<size_t, double>& l, const std::pair<size_t, double>& r) {
            return l.second > r.second; });
    
    size_t begin = std::min(gaps[0].first, gaps[1].first);
    size_t before_end = std::max(gaps[0].first, gaps[1].first);
    ++begin;
    std::copy(planes_its.begin() + begin, planes_its.begin() + before_end + 1, std::back_inserter(rundist_planes_its));


    double z_start = rundist_planes_its.begin()->first.c();
    double z_end = std::prev(rundist_planes_its.end())->first.c();
    for(auto& el: planes_its) {
        if (el.first.c() > z_start) {
            up_rundist_planes_its.push_back(el);
        } else if (el.first.c() < z_end) {
            low_rundist_planes_its.push_back(el);
        }
    }

    return;
}

double calc_height(std::vector<std::pair<Plane, Face_iterator>>& v1) {
    double min_z = std::numeric_limits<double>::max();
    double max_z = std::numeric_limits<double>::min();
    for(auto& el: v1) {
        auto begin = el.second->facet_begin();
        size_t sz = el.second->size();
        for(size_t i = 0; i < sz; ++i, ++begin) {
            double z = begin->vertex()->point().z();
            if (z < min_z) {
                min_z = z;
            }
            if (z > max_z) {
                max_z = z;
            }
        }
    }
    return max_z - min_z;
}

inline Plane unit_plane_equation(Point_3 p1, Point_3 p2, Point_3 p3) {
    Plane plane(p1, p2, p3);
    double norm = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
    Plane norm_plane(plane.a()/norm, plane.b()/norm, plane.c()/norm, plane.d()/norm);
    return norm_plane;
}


void get_plane_equations_from_polyhedron(Polyhedron& input_poly
            , std::vector<std::pair<Plane, Face_iterator>>& planes_its) {
    auto end = input_poly.facets_end();
    for(auto it = input_poly.facets_begin(); it != end; ++it) {
        auto h = it->halfedge();
        Point_3 p1 = h->vertex()->point();
        Point_3 p2 = h->next()->vertex()->point();
        Point_3 p3 = h->next()->next()->vertex()->point();
        while (CGAL::collinear(p1, p2, p3)) {
            h = h->next();
            p1 = h->vertex()->point();
            p2 = h->next()->vertex()->point();
            p3 = h->next()->next()->vertex()->point();
        }

        Plane plane = unit_plane_equation(p1, p2, p3);
        planes_its.push_back({plane, it});
    }
}

#define _DTOR_   0.0174532925199432957692         //  PI / 180
inline double Rad2Deg(double a) { return a / _DTOR_; }
inline double Deg2Rad(double a) { return a * _DTOR_; }
 
double Azimuth(const Point_3& v) {
    double a = atan2(v.x(), v.y());
    return (a >= 0) ? a : a + 2*M_PI;
}

double Slope(const Point_3 &v) {
    double d = sqrt(v.x()*v.x() + v.y()*v.y());
    return atan2(v.z(), d);
}

Point_3 VectorFrom(double azimuth, double slope, double length) {
    double c = cos(slope);
    Point_3 v(length*sin(azimuth)*c, length*cos(azimuth)*c, length*sin(slope));
    return v;
}

void transform_height_convex(std::vector<Plane>& planes, std::vector<std::pair<Plane, Face_iterator>>& to_change
        , std::unordered_set<Vertex_handle>& rundist_vertices, double delta, bool pavilion) {
    std::vector<Point_3> point_cloud;
    for(auto& el: to_change) {
        auto begin = std::get<1>(el)->facet_begin();
        auto end = rundist_vertices.end();
        size_t sz = std::get<1>(el)->size();
        for(size_t i = 0; i < sz; ++i, ++begin) {
            auto vh = begin->vertex();
            Point_3 p = vh->point();
            if (auto it = rundist_vertices.find(vh); it == end) {
                point_cloud.push_back(Point_3(p.x(), p.y(), p.z() + delta));
            } else {
                point_cloud.push_back(p);
            }
        }
    }
    Polyhedron upper_polly;
    CGAL::convex_hull_3(point_cloud.begin(), point_cloud.end(), upper_polly);
    std::vector<std::pair<Vector_3, Face_iterator>> upper_poly_faces;
    auto end = upper_polly.facets_end();
    for(auto it = upper_polly.facets_begin(); it != end; ++it) {
        auto hv = it->halfedge();
        Point_3 p1 = hv->vertex()->point();
        Point_3 p2 = hv->next()->vertex()->point();
        Point_3 p3 = hv->next()->next()->vertex()->point();
        Vector_3 unit_v = CGAL::unit_normal(p1, p2, p3); 
        upper_poly_faces.push_back({unit_v, it});
    }
    std::sort(upper_poly_faces.begin(), upper_poly_faces.end(), [](const std::pair<Vector_3, Face_iterator>& lf
                , const std::pair<Vector_3, Face_iterator>& rf) { return lf.first.z() > rf.first.z(); });
    size_t sz = upper_poly_faces.size();
    size_t index = 0;
    for(size_t i = 0; i < sz - 1; ++i) {
        if (upper_poly_faces[i].first.z() < 0) {
            index = i;
            break;
        }
    }
    size_t begin = (pavilion == true) ? 0 : index;
    size_t stop = (pavilion == true) ? index : sz;
    for(size_t i = begin; i < stop; ++i) {
        auto hv = upper_poly_faces[i].second->halfedge();
        Point_3 p1 = hv->vertex()->point();
        Point_3 p2 = hv->next()->vertex()->point();
        Point_3 p3 = hv->next()->next()->vertex()->point();
        Plane plane = unit_plane_equation(p1, p2, p3);
        planes.push_back(plane);
    }
}

void apply_changes(std::vector<std::pair<Plane, Face_iterator>>& planes, double parametr_of_azimuth
        , double parametr_of_slope, double parametr_of_shift, double sigma_shift
        , double sigma_slope, double sigma_azimuth, bool is_change_azimuth, bool is_change_slope
        , bool is_change_shift, std::string distribution) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_uni_azimuth(-parametr_of_azimuth, parametr_of_azimuth);
    std::uniform_real_distribution<> dis_uni_slope(-parametr_of_slope, parametr_of_slope);
    std::uniform_real_distribution<> dis_uni_shift(-parametr_of_shift, parametr_of_shift);
    std::normal_distribution<> dis_normal_shift(0, sigma_shift);
    std::normal_distribution<> dis_normal_slope(0, sigma_slope);
    std::normal_distribution<> dis_normal_azimuth(0, sigma_azimuth);

    srand(time(NULL));
    for(auto& el: planes) {
        Plane plane = std::get<0>(el);
        Point_3 vec(plane.a(), plane.b(), plane.c());
        double d = plane.d();
        double azimuth = Azimuth(vec);
        double slope = Slope(vec);
        if (is_change_azimuth == true) {
            if (distribution == "uniform") {
                azimuth += dis_uni_azimuth(gen);
            } else {
                double l = dis_normal_azimuth(gen);
                if (std::fabs(l) > 2 * sigma_azimuth) {
                    l = (l > 0) ? 2 * sigma_azimuth : -2 * sigma_azimuth;
                }
                azimuth += l;
            }   
        }
        if (is_change_slope == true) {
            if (distribution == "uniform") {
                slope += dis_uni_slope(gen);
            } else {
                double l = dis_normal_slope(gen);
                if (std::fabs(l) > 2 * sigma_slope) {
                    l = (l > 0) ? 2 * sigma_slope : -2 * sigma_slope;
                }
                slope += l;
            }   
        }
        vec = VectorFrom(azimuth, slope, 1);
        if (is_change_shift == true) {
            if (distribution == "uniform") {
                d += dis_uni_shift(gen);
            } else {
                double l = dis_normal_shift(gen);
                if (std::fabs(l) > 2 * sigma_shift) {
                    l = (l > 0) ? 2 * sigma_shift : -2 * sigma_shift;
                }
                d += l;
            }
        }
        Plane plane_new(vec.x(), vec.y(), vec.z(), d);
        std::get<0>(el) = plane_new;
    }
}


int main(int argc, char* argv[]) {
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    std::string input_file;
    std::string type_of_input_file;
    std::string output_file;
    std::string type_of_output_file;
    bool is_change_shift;
    bool is_change_azimuth;
    bool is_change_slope;
    bool is_change_height_crown;
    bool is_change_height_pavilion;
    double parametr_of_shift;
    double parametr_of_azimuth;
    double parametr_of_slope;
    double parametr_of_height_crown;
    double parametr_of_height_pavilion;
    bool is_change_crown;
    bool is_change_pavilion;
    bool is_change_rundist;
    std::string distribution;
    double sigma_shift;
    double sigma_azimuth;
    double sigma_slope;
    boost::program_options::options_description desc("All options");
    desc.add_options()
        ("input_file", boost::program_options::value<std::string>(&input_file) 
         -> default_value("model_in.ply"), "file with initial model")
        ("type_of_input_file", boost::program_options::value<std::string>(&type_of_input_file) 
         -> default_value("ply"), "txt, ply or obj")
        ("output_file", boost::program_options::value<std::string>(&output_file)
         -> default_value("model_out.ply"), "file to save model in specila format")
        ("type_of_output_file", boost::program_options::value<std::string>(&type_of_output_file) 
         -> default_value("ply"), "ply or obj")

        ("is_change_shift", boost::program_options::value<bool>(&is_change_shift)
         -> default_value(false), "do you want to shift a parametr d")
        ("is_change_azimuth", boost::program_options::value<bool>(&is_change_azimuth)
         -> default_value(false), "do you want to change azimuth")
        ("is_change_slope", boost::program_options::value<bool>(&is_change_slope)
         -> default_value(false), "do you want to change slope")
        ("is_change_height_crown", boost::program_options::value<bool>(&is_change_height_crown)
         -> default_value(false), "do you want to change height of crown")
        ("is_change_height_pavilion", boost::program_options::value<bool>(&is_change_height_pavilion)
         -> default_value(false), "do you want to change height of pavilion")

        ("parametr_of_shift", boost::program_options::value<double>(&parametr_of_shift)
         -> default_value(0.0001), " d += parametr_of_shift")
        ("parametr_of_azimuth", boost::program_options::value<double>(&parametr_of_azimuth)
         -> default_value(1), " azimuth += parametr_of_azimuth (in degrees)")
        ("parametr_of_slope", boost::program_options::value<double>(&parametr_of_slope)
         -> default_value(1), " slope += parametr_of_slope (in degrees)")
        ("parametr_of_height_crown", boost::program_options::value<double>(&parametr_of_height_crown)
         -> default_value(105), " percents of height of an old height")
        ("parametr_of_height_pavilion", boost::program_options::value<double>(&parametr_of_height_pavilion)
         -> default_value(105), " percents of height of an old height")

        ("is_change_crown", boost::program_options::value<bool>(&is_change_crown)
         -> default_value(false), "do you want to apply changes to crown")
        ("is_change_pavilion", boost::program_options::value<bool>(&is_change_pavilion)
         -> default_value(false), "do you want to apply changes to pavilion")
        ("is_change_rundist", boost::program_options::value<bool>(&is_change_rundist)
         -> default_value(false), "do you want to apply changes to rundist")

        ("distribution", boost::program_options::value<std::string>(&distribution)
         -> default_value("uniform"), "uniform or normal")
        ("sigma_shift", boost::program_options::value<double>(&sigma_shift)
         -> default_value(0.0001), "dispersion in normal distribution for shift")
        ("sigma_azimuth", boost::program_options::value<double>(&sigma_azimuth)
         -> default_value(1), "dispersion in normal distribution for azimuth (in derees)")
        ("sigma_slope", boost::program_options::value<double>(&sigma_slope)
         -> default_value(1), "dispersion in normal distribution for slope (in degrees)")

        ("help", "produce help message")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    sigma_azimuth = Deg2Rad(sigma_azimuth);
    sigma_slope = Deg2Rad(sigma_slope);
    parametr_of_azimuth = Deg2Rad(parametr_of_azimuth);
    parametr_of_slope = Deg2Rad(parametr_of_slope);

    Polyhedron input_poly;
    std::vector<std::pair<Plane, Face_iterator>> planes_its;

    bool succes = false;
    if (type_of_input_file == "ply") {
        succes = CGAL::IO::read_PLY(input_file, input_poly);
        if (succes == false) {
            std::cerr << "Cannot read ply polyhedron from: " << input_file << std::endl;
            return -1;
        }
        get_plane_equations_from_polyhedron(input_poly, planes_its);

    } else if (type_of_input_file == "obj") {
        succes = CGAL::IO::read_OBJ(input_file, input_poly);
        if (succes == false) {
            std::cerr << "Cannot read obj polyhedron from: " << input_file << std::endl;
            return -1;
        }
        get_plane_equations_from_polyhedron(input_poly, planes_its);

    } else if (type_of_input_file == "txt") {
        std::vector<Plane> planes;
        succes = read_planes_from_file(planes, input_file);
        if (succes == false) {
            std::cerr << "Cannot read planes from txt file from: " << input_file << std::endl;
            return -1;
        }
        CGAL::halfspace_intersection_with_constructions_3(planes.begin(), planes.end(), input_poly); 
        get_plane_equations_from_polyhedron(input_poly, planes_its);
    } else {
        std::cerr << "Wrong type of file" << std::endl;
        return -1;
    }
    std::vector<std::pair<Plane, Face_iterator>> rundist_planes_its;
    std::vector<std::pair<Plane, Face_iterator>> up_rundist_planes_its;
    std::vector<std::pair<Plane, Face_iterator>> low_rundist_planes_its;

    find_rundist_up_low(planes_its, rundist_planes_its, up_rundist_planes_its, low_rundist_planes_its);

#ifdef DEBUG
    print_polyhedron_ply(input_poly, "debug_initial_polyhedron.ply");
    print_faces_as_polyhedron_ply(rundist_planes_its, 0, rundist_planes_its.size(), "debug_rundist.ply");
    print_faces_as_polyhedron_ply(up_rundist_planes_its, 0, up_rundist_planes_its.size(), "debug_up_rundist.ply");
    print_faces_as_polyhedron_ply(low_rundist_planes_its, 0, low_rundist_planes_its.size(), "debug_low_rundist.ply");
#endif


    std::unordered_set<Vertex_handle> rundist_vertices;
    for(auto& el: rundist_planes_its) {
        auto begin = el.second->facet_begin();
        size_t sz = el.second->size();
        for(size_t i = 0; i < sz; ++i, ++begin) {
            rundist_vertices.insert(begin->vertex());
        }
    }

   
    std::vector<Plane> planes;

    if (is_change_height_pavilion == true) {
        double height_pavilion = calc_height(up_rundist_planes_its);
        double delta = height_pavilion * (parametr_of_height_pavilion/100 - 1);
        transform_height_convex(planes, up_rundist_planes_its, rundist_vertices, delta, true);
    }

    if (is_change_height_crown == true) {
        double height_crown = calc_height(low_rundist_planes_its);
        double delta = height_crown * (parametr_of_height_crown/100 - 1);
        transform_height_convex(planes, low_rundist_planes_its, rundist_vertices, -delta, false);
    }


    if(is_change_pavilion == true) {
        apply_changes(up_rundist_planes_its, parametr_of_azimuth, parametr_of_slope, parametr_of_shift, sigma_shift
                , sigma_slope, sigma_azimuth, is_change_azimuth, is_change_slope, is_change_shift, distribution);
    }

    if(is_change_crown == true) {
        apply_changes(low_rundist_planes_its, parametr_of_azimuth, parametr_of_slope, parametr_of_shift, sigma_shift
                , sigma_slope, sigma_azimuth, is_change_azimuth, is_change_slope, is_change_shift, distribution);
    }

    if(is_change_rundist == true) {
        apply_changes(rundist_planes_its, parametr_of_azimuth, parametr_of_slope, parametr_of_shift, sigma_shift
                , sigma_slope, sigma_azimuth, is_change_azimuth, is_change_slope, is_change_shift, distribution);
    }

    if (is_change_height_pavilion == false) {
        for(auto& el: up_rundist_planes_its) {
            planes.push_back(std::get<0>(el));
        }
    }
    for(auto& el: rundist_planes_its) {
        planes.push_back(el.first);
    }
    if (is_change_height_crown == false) {
        for(auto& el: low_rundist_planes_its) {
            planes.push_back(std::get<0>(el));
        }
    }
   
    Polyhedron output_poly;
    CGAL::halfspace_intersection_with_constructions_3(planes.begin(), planes.end(), output_poly); 

    if (type_of_output_file == "ply") {
        std::ofstream of1(output_file);
        of1 << std::fixed << std::setprecision(12);
        CGAL::IO::write_PLY(of1, output_poly);
    } else if (type_of_output_file == "obj") {
        std::ofstream of2(output_file);
        of2 << std::fixed << std::setprecision(12);
        CGAL::IO::write_OBJ(of2, output_poly);
    }

    return 0;
}




    

    

