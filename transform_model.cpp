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
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;
typedef K::Point_3                                          Point_3;
typedef K::Plane_3                                          Plane;
typedef K::Vector_3                                         Vector_3;
typedef Polyhedron::Face_iterator                           Face_iterator; 
typedef Polyhedron::Face_handle                             Face_handle;
typedef Polyhedron::Vertex_handle                           Vertex_handle;



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

void find_rundist_and_other_parts(std::vector<std::pair<Plane, Face_iterator>>& planes_its
        , std::vector<std::pair<Plane, Face_iterator>>& rundist_planes_its
        , std::vector<std::tuple<Plane, Face_iterator, Point_3, Point_3>>& up_rundist_planes_its
        , std::vector<std::pair<Plane, Face_iterator>>& up_up_rundist_planes_its
        , std::vector<std::tuple<Plane, Face_iterator, Point_3, Point_3>>& low_rundist_planes_its
        , std::vector<std::pair<Plane, Face_iterator>>& low_low_rundist_planes_its
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

    std::map<Face_handle, std::tuple<Plane, Face_iterator, Point_3, Point_3>> up_rundist;
    std::map<Face_handle, std::tuple<Plane, Face_iterator, Point_3, Point_3>> low_rundist;


    double z_start = rundist_planes_its.begin()->first.c();
    double z_end = std::prev(rundist_planes_its.end())->first.c();
    for(auto& el: rundist_planes_its) {
        auto begin = el.second->facet_begin();
        size_t sz = el.second->size(); 
        for(size_t j = 0; j < sz; ++j, ++begin) {
            Point_3 p1 = begin->vertex()->point();
            Point_3 p2 = begin->opposite()->vertex()->point();
            auto op_face = begin->opposite()->facet();
            auto h = op_face->halfedge();
            Plane plane = Plane(h->vertex()->point(), h->next()->vertex()->point()
                    , h->next()->next()->vertex()->point());
            double norm = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
            Plane norm_plane(plane.a()/norm, plane.b()/norm, plane.c()/norm, plane.d()/norm);
            std::tuple<Plane, Face_iterator, Point_3, Point_3> p(norm_plane, op_face, p1, p2);
            if (norm_plane.c() > z_start) {
                up_rundist.insert(std::pair(op_face, p));
            } 
            if (norm_plane.c() < z_end) {
                low_rundist.insert(std::pair(op_face, p));
            }
        }
    }

    auto end = planes_its.begin() + begin + 1;
    for(auto begin = planes_its.begin(); begin != end; ++begin) {
        if(auto it = up_rundist.find(begin->second); it == up_rundist.end()) {
            up_up_rundist_planes_its.push_back(*begin);
        }
    }
    end = planes_its.end();
    for(auto begin = planes_its.begin() + before_end + 1; begin != end; ++begin) {
        if(auto it = low_rundist.find(begin->second); it == low_rundist.end()) {
            low_low_rundist_planes_its.push_back(*begin);
        }
    }

    for(auto& el: up_rundist) {
        up_rundist_planes_its.push_back(el.second);
    }
    for(auto& el: low_rundist) {
        low_rundist_planes_its.push_back(el.second);
    }

    return;
}

double calc_height(std::vector<std::pair<Plane, Face_iterator>>& v1
        , std::vector<std::tuple<Plane, Face_iterator, Point_3, Point_3>>& v2) {
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
    for(auto& el: v2) {
        auto begin = std::get<1>(el)->facet_begin();
        size_t sz = std::get<1>(el)->size();
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

void get_plane_equations_from_polyhedron(Polyhedron& input_poly
            , std::vector<std::pair<Plane, Face_iterator>>& planes_its) {
    auto end = input_poly.facets_end();
    for(auto it = input_poly.facets_begin(); it != end; ++it) {
        auto h = it->halfedge();
        Point_3 p1 = h->vertex()->point();
        Point_3 p2 = h->next()->vertex()->point();
        Point_3 p3 = h->next()->next()->vertex()->point();
        Plane plane(p1, p2, p3);
        double norm = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
        Plane norm_plane(plane.a()/norm, plane.b()/norm, plane.c()/norm, plane.d()/norm);
        planes_its.push_back({norm_plane, it});
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




int main(int argc, char* argv[]) {
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
    std::vector<std::tuple<Plane, Face_iterator, Point_3, Point_3>> up_rundist_palnes_its;
    std::vector<std::pair<Plane, Face_iterator>> up_up_rundist_palnes_its;
    std::vector<std::tuple<Plane, Face_iterator, Point_3, Point_3>> low_rundist_palnes_its;
    std::vector<std::pair<Plane, Face_iterator>> low_low_rundist_palnes_its;

    find_rundist_and_other_parts(planes_its, rundist_planes_its, up_rundist_palnes_its
            , up_up_rundist_palnes_its, low_rundist_palnes_its, low_low_rundist_palnes_its);


    std::unordered_set<Vertex_handle> rundist_vertices;
    for(auto& el: rundist_planes_its) {
        auto begin = el.second->facet_begin();
        size_t sz = el.second->size();
        for(size_t i = 0; i < sz; ++i, ++begin) {
            rundist_vertices.insert(begin->vertex());
        }
    }

   
    double height_pavilion = 0, height_crown =0;

    height_pavilion = calc_height(up_up_rundist_palnes_its, up_rundist_palnes_its);
    height_crown = calc_height(low_low_rundist_palnes_its, low_rundist_palnes_its);

    std::vector<Plane> planes;

    if (is_change_height_pavilion == true) {
        double delta = height_pavilion * (parametr_of_height_pavilion/100 - 1);
        for(auto& el: up_rundist_palnes_its) {
            std::vector<Point_3> part_of_rundist;
            std::vector<Point_3> not_rundist;
            auto begin = std::get<1>(el)->facet_begin();
            auto vh = begin->vertex();
            if (auto it = rundist_vertices.find(vh); it != rundist_vertices.end()) {
                auto end = rundist_vertices.end();
                while(it != end) {
                    ++begin;
                    vh = begin->vertex();
                    it = rundist_vertices.find(vh);
                }
            }
            size_t sz = std::get<1>(el)->size();
            for(size_t i = 0; i < sz; ++i, ++begin) {
                auto vh = begin->vertex();
                if (auto it = rundist_vertices.find(vh); it != rundist_vertices.end()) {
                    part_of_rundist.push_back(vh->point());
                } else {
                    not_rundist.push_back(vh->point());
                }
            }

            size_t sz1 = part_of_rundist.size(), sz2 = not_rundist.size();
            for(size_t i = 0; i < sz1 - 1; ++i) {
                Point_3 p1 = part_of_rundist[i];
                Point_3 p2 = part_of_rundist[i + 1];
                for(size_t j = 0; j < sz2; ++j) {
                    Point_3 p3 = not_rundist[j];
                    p3 = Point_3(p3.x(), p3.y(), p3.z() + delta);
                    Plane plane(p3, p1, p2);
                    double norm = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
                    Plane norm_plane(plane.a()/norm, plane.b()/norm, plane.c()/norm, plane.d()/norm);
                    planes.push_back(norm_plane);
                }
            }
        }
        for(auto& el: up_up_rundist_palnes_its) {
            auto h = el.second->halfedge();
            Point_3 p1 = h->vertex()->point();
            Point_3 p2 = h->next()->vertex()->point();
            Point_3 p3 = h->next()->next()->vertex()->point();
            p1 = Point_3(p1.x(), p1.y(), p1.z() + delta);
            p2 = Point_3(p2.x(), p2.y(), p2.z() + delta);
            p3 = Point_3(p3.x(), p3.y(), p3.z() + delta);
            Plane plane(p1, p2, p3);
            double norm = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
            Plane norm_plane(plane.a()/norm, plane.b()/norm, plane.c()/norm, plane.d()/norm);
            el.first = norm_plane;
        }
    }

    if (is_change_height_crown == true) {
        double delta = height_crown * (parametr_of_height_crown/100 - 1);
        std::cout << delta << std::endl;
        for(auto& el: low_rundist_palnes_its) {
            std::vector<Point_3> part_of_rundist;
            std::vector<Point_3> not_rundist;
            auto begin = std::get<1>(el)->facet_begin();
            auto vh = begin->vertex();
            if (auto it = rundist_vertices.find(vh); it != rundist_vertices.end()) {
                auto end = rundist_vertices.end();
                while(it != end) {
                    ++begin;
                    vh = begin->vertex();
                    it = rundist_vertices.find(vh);
                }
            }
            size_t sz = std::get<1>(el)->size();
            for(size_t i = 0; i < sz; ++i, ++begin) {
                auto vh = begin->vertex();
                if (auto it = rundist_vertices.find(vh); it != rundist_vertices.end()) {
                    part_of_rundist.push_back(vh->point());
                } else {
                    not_rundist.push_back(vh->point());
                }
            }

            size_t sz1 = part_of_rundist.size(), sz2 = not_rundist.size();
            for(size_t i = 0; i < sz1 - 1; ++i) {
                Point_3 p1 = part_of_rundist[i];
                Point_3 p2 = part_of_rundist[i + 1];
                for(size_t j = 0; j < sz2; ++j) {
                    Point_3 p3 = not_rundist[j];
                    p3 = Point_3(p3.x(), p3.y(), p3.z() - delta);
                    std::cout << p3.z() << std::endl;
                    Plane plane(p3, p1, p2);
                    double norm = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
                    Plane norm_plane(plane.a()/norm, plane.b()/norm, plane.c()/norm, plane.d()/norm);
                    planes.push_back(norm_plane);
                }
            }
        }
        for(auto& el: low_low_rundist_palnes_its) {
            auto h = el.second->halfedge();
            Point_3 p1 = h->vertex()->point();
            Point_3 p2 = h->next()->vertex()->point();
            Point_3 p3 = h->next()->next()->vertex()->point();
            p1 = Point_3(p1.x(), p1.y(), p1.z() - delta);
            p2 = Point_3(p2.x(), p2.y(), p2.z() - delta);
            p3 = Point_3(p3.x(), p3.y(), p3.z() - delta);
            Plane plane(p1, p2, p3);
            double norm = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
            Plane norm_plane(plane.a()/norm, plane.b()/norm, plane.c()/norm, plane.d()/norm);
            el.first = norm_plane;
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_uni_azimuth(-parametr_of_azimuth, parametr_of_azimuth);
    std::uniform_real_distribution<> dis_uni_slope(-parametr_of_slope, parametr_of_slope);
    std::uniform_real_distribution<> dis_uni_shift(-parametr_of_shift, parametr_of_shift);
    std::normal_distribution<> dis_normal_shift(0, sigma_shift);
    std::normal_distribution<> dis_normal_slope(0, sigma_slope);
    std::normal_distribution<> dis_normal_azimuth(0, sigma_azimuth);


    srand(time(NULL));
    if(is_change_pavilion == true) {
        for(auto& el: up_rundist_palnes_its) {
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
        for(auto& el: up_up_rundist_palnes_its) {
            Plane plane = el.first;
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
            el.first = plane_new;
        }
    }

    if(is_change_crown == true) {
        for(auto& el: low_rundist_palnes_its) {
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
        for(auto& el: low_low_rundist_palnes_its) {
            Plane plane = el.first;
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
            el.first = plane_new;
        }
    }

    if(is_change_rundist == true) {
        for(auto& el: rundist_planes_its) {
            Plane plane = el.first;
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
            el.first = plane_new;
        }
    }

    for(auto& el: up_up_rundist_palnes_its) {
        planes.push_back(el.first);
    }
    if (is_change_height_pavilion == false) {
        for(auto& el: up_rundist_palnes_its) {
            planes.push_back(std::get<0>(el));
        }
    }
    for(auto& el: rundist_planes_its) {
        planes.push_back(el.first);
    }
    if (is_change_height_crown == false) {
        for(auto& el: low_rundist_palnes_its) {
            planes.push_back(std::get<0>(el));
        }
    }
    for(auto& el: low_low_rundist_palnes_its) {
        planes.push_back(el.first);
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




    

    

