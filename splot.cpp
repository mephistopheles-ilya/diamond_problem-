#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Polyhedron_3_fwd.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/IO/STL.h>
#include <CGAL/IO/PLY.h>

#include <boost/program_options.hpp>

#include <cmath>
#include <cstdio>
#include <list>
#include <string>
#include <istream>
#include <random>

#include "include/point3d.hpp"

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Polyhedron_3<K>                                 Polyhedron_3;
typedef K::Plane_3                                            Plane;
typedef K::Point_3                                            Point;
typedef CGAL::Surface_mesh<Point>                             Surface_mesh;

#define _DTOR_   0.0174532925199432957692         //  PI / 180
inline double Rad2Deg(double a) { return a / _DTOR_; }
inline double Deg2Rad(double a) { return a * _DTOR_; }
 
double Azimuth(const Point3D& v) {
    double a = atan2(v.x, v.y);
    return (a >= 0) ? a : a + 2*M_PI;
}

double Slope(const Point3D &v) {
    double d = sqrt(v.x*v.x + v.y*v.y);
    return atan2(v.z, d);
}

Point3D VectorFrom(double azimuth, double slope, double length) {
    double c = cos(slope);
    Point3D v(length*sin(azimuth)*c, length*cos(azimuth)*c, length*sin(slope));
    return v;
}

int sign(double number) {
    if (number > 0) {
        return 1;
    }
    if (number < 0) {
        return -1;
    }
    return 0;
}



void read_planes(std::istream& is, std::list<Plane>& planes, bool is_change_shift,  double parametr_of_shift
        , bool is_change_azimuth, double parametr_of_azimuth, bool is_change_slope, double parametr_of_slope
        , int sign_of_c, std::string distribution, double sigma) {
    std::string line;
    //std::string key("# plane coeff");
    std::string stop("#");
    std::string key("indices_of_vertices");
    for(;;) {
        std::getline(is, line);
        if (line.find(key) != std::string::npos) 
            break;
    }
    std::getline(is, line);
    int counter = 1;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_uni_azimuth(-parametr_of_azimuth, parametr_of_azimuth);
    std::uniform_real_distribution<> dis_uni_slope(-parametr_of_slope, parametr_of_slope);
    std::uniform_real_distribution<> dis_uni_shift(-parametr_of_shift, parametr_of_shift);
    std::normal_distribution<> dis_normal_azimuth(0, sigma);
    std::normal_distribution<> dis_normal_slope(0, sigma);
    std::normal_distribution<> dis_normal_shift(0, sigma);

    srand(time(NULL));
    while(line.find(stop) == std::string::npos) {
        double a, b, c, d;
        int unused1, unused2;
        sscanf(line.c_str(), "%d %d %lf %lf %lf %lf", &unused1, &unused2, &a, &b, &c, &d);
#if 0
        if (c > 0) {
            typename K::Plane_3 plane(a, b , c, d + 0.001);
            planes.push_back(plane);
        } else {
            typename K::Plane_3 plane(a, b, c, d);
            planes.push_back(plane);
        }
#endif
#if 0
       if (sign_of_c == sign(c)) {
            Point3D vec(a, b, c);
            double azimuth = Azimuth(vec);
            double slope = Slope(vec);
            if (is_change_azimuth == true) {
                azimuth += parametr_of_azimuth;
            }
            if (is_change_slope == true) {
                slope += parametr_of_slope;
            }
            vec = VectorFrom(azimuth, slope, 1);
            if (is_change_shift == true) {
                typename K::Plane_3 plane(vec.x, vec.y, vec.z, d + parametr_of_shift);
                planes.push_back(plane);
            } else {
                typename K::Plane_3 plane(vec.x, vec.y, vec.z, d);
                planes.push_back(plane);
            }
        } else {
            typename K::Plane_3 plane(a, b, c, d);
            planes.push_back(plane);
        }
#endif
        if (sign_of_c == sign(c)) {
            Point3D vec(a, b, c);
            double azimuth = Azimuth(vec);
            double slope = Slope(vec);
            if (is_change_azimuth == true) {
                if (distribution == "uniform") {
                    azimuth += dis_uni_azimuth(gen);
                } else {
                    double l = dis_normal_azimuth(gen);
                    if (std::fabs(l) > parametr_of_azimuth) {
                        l = (l > 0) ? parametr_of_azimuth : -parametr_of_azimuth;
                    }
                    azimuth += l;
                }   
            }
            if (is_change_slope == true) {
                if (distribution == "uniform") {
                    slope += dis_uni_slope(gen);
                } else {
                    double l = dis_normal_slope(gen);
                    if (std::fabs(l) > parametr_of_slope) {
                        l = (l > 0) ? parametr_of_slope : -parametr_of_slope;
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
                    if (std::fabs(l) > parametr_of_shift) {
                        l = (l > 0) ? parametr_of_shift : -parametr_of_shift;
                    }
                    d += l;
                }   
                typename K::Plane_3 plane(vec.x, vec.y, vec.z, d);
                planes.push_back(plane);
            } else {
                typename K::Plane_3 plane(vec.x, vec.y, vec.z, d);
                planes.push_back(plane);
            }
        } else {
            typename K::Plane_3 plane(a, b, c, d);
            planes.push_back(plane);
        }

        std::getline(is, line);
        std::getline(is, line);
        ++counter;
    }
}
 
int main (int argc, char* argv[]) {
    std::string initial_model;
    std::string outputfile;
    bool is_change_shift;
    bool is_change_azimuth;
    bool is_change_slope;
    double parametr_of_shift;
    double parametr_of_azimuth;
    double parametr_of_slope;
    int sign_of_c;
    std::string distribution;
    double sigma;
    boost::program_options::options_description desc("All options");
    desc.add_options()
        ("initial_model", boost::program_options::value<std::string>(&initial_model) 
         -> default_value("InitialModels/InitialModel5"), "file with initial model")
        ("output", boost::program_options::value<std::string>(&outputfile)
         -> default_value("tmp.ply"), "file to save model in ply format")
        ("is_change_shift", boost::program_options::value<bool>(&is_change_shift)
         -> default_value(false), "do you want to shift a parametr d")
        ("is_change_azimuth", boost::program_options::value<bool>(&is_change_azimuth)
         -> default_value(false), "do you want to change azimuth")
        ("is_change_slope", boost::program_options::value<bool>(&is_change_slope)
         -> default_value(false), "do you want to change slope")
        ("parametr_of_shift", boost::program_options::value<double>(&parametr_of_shift)
         -> default_value(0.0001), " d += parametr_of_shift")
        ("parametr_of_azimuth", boost::program_options::value<double>(&parametr_of_azimuth)
         -> default_value(0.0001), " azimuth += parametr_of_azimuth")
        ("parametr_of_slope", boost::program_options::value<double>(&parametr_of_slope)
         -> default_value(0.0001), " slope += parametr_of_slope")
        ("sign_of_c", boost::program_options::value<int>(&sign_of_c)
         -> default_value(1), "sign of coefficient c of palnes which will be changed")
        ("distribution", boost::program_options::value<std::string>(&distribution)
         -> default_value("uniform"), "distribution in changind vectors of facets")
        ("sigma", boost::program_options::value<double>(&sigma)
         -> default_value(0.0001), "dispearsion in normal distribution")

        ("help", "produce help message")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    std::list<Plane> planes;
    std::ifstream is(initial_model);
    read_planes(is, planes , is_change_shift, parametr_of_shift
            , is_change_azimuth, parametr_of_azimuth, is_change_slope, parametr_of_slope
            , sign_of_c, distribution, sigma); 

 
    Polyhedron_3 chull;

    CGAL::halfspace_intersection_with_constructions_3(planes.begin(), planes.end(), chull);

    //std::ofstream of1("tmp.stl");
    std::ofstream of2(outputfile);
    //CGAL::IO::write_STL(of1, chull);
    CGAL::IO::write_PLY(of2, chull);

#if 0
    for(auto it = chull.facets_begin(); it != chull.facets_end(); ++it) {
        auto v = it->facet_begin();
        for(int i = 0; i < it->size(); ++i) {
            std::cout << v->vertex()->point() << std::endl;
            ++v;
        }
        std::cout << it->facet_begin()->vertex()->point() << std::endl;
        std::cout << std::endl << std::endl;
    }
#endif

    return 0;
}
