#include "include/point2d.h"
#include "include/point3d.h"
#include "include/edge.h"
#include "include/rw_functions.h"
#include "include/geom_functions.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numbers>
#include <cmath>

#include <boost/program_options.hpp>

int main(int argc, char* argv[]) {
//PARSING OF COMMAND LINE
    int projections;
    std::string file;
    std::string projection_method;
    boost::program_options::options_description desc("All options");
    desc.add_options()
        ("projections", boost::program_options::value<int>(&projections) -> default_value(30)
         , "how many projections do you want to create")
        ("init_file", boost::program_options::value<std::string>(&file) 
         -> default_value("InitialModel"), "file with initial model")
        ("pr_method", boost::program_options::value<std::string>(&projection_method)
         -> default_value("hull"), "method to make projection : hull, obvious, graph")
        ("help", "produce help message")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    std::cout << "projections = " << projections << '\n';
    std::cout << "file = " << file << '\n';
    std::cout << "method = " << projection_method << std::endl;


    std::ifstream in(file);
    if(in.is_open()) {
        std::cout << "Input files are opened" << std::endl;
    } else {
        std::cout << "Wrong files" << std::endl;
        throw std::runtime_error("Cannot open or read file\n");
    }
    in >> std::fixed >> std::setprecision(PRECISION);

//READING STRUCTURES FROM FILE
    std::vector<Point3D> arr_points3d;
    std::vector<Point3D> arr_norm3d;
    std::vector<Edge> arr_edges3d;
    auto [num_vertices, num_facest, num_edges] = read_structures_from_file(in
            , arr_points3d, arr_norm3d, arr_edges3d); 
    std::cout << "num_vertices = " << num_vertices << ' ' << arr_points3d.size() << '\n';
    std::cout << "num_facest = " << num_facest << ' ' << arr_norm3d.size() << '\n';
    std::cout << "num_edges = " << num_edges << ' ' << arr_edges3d.size() <<  std::endl;

//CREATING PROJECTIONS
    for (int step = 0; step <= projections; ++step) { 
        boost::geometry::model::polygon<Point2D, false, true, std::vector> hull;
        double angle = step * std::numbers::pi_v<double>/projections;
//CHOOSING METHOD TO CREATE PRIJECTIONS
        if (projection_method == "hull") {
            projection_convex_hull(arr_points3d, hull, angle);
        } else if (projection_method == "obvious") {
            projection_square_complexity(arr_points3d, arr_norm3d, arr_edges3d, hull, angle);
        } else if (projection_method == "graph") { 
            projection_graph(arr_points3d, arr_norm3d, arr_edges3d, hull, angle);
        }

//CHECK IF PROJECTION IS VALID
        boost::geometry::correct(hull);
        std::cout << "Hull size : " << hull.outer().size() << std::endl;
        std::string message;
        if(!boost::geometry::is_valid(hull, message)) {
            std::cerr << "convex hull is Not valid " << message << std::endl;
        }

//TRYING TO DELETE EXTRA PART IN THE BOTTOM
        Point2D p_min_y = find_minimal_y(hull);

        double min_y_coord = p_min_y.y;
        auto [min_x_coord, max_x_coord] = find_borders_for_x(p_min_y, hull);

        boost::geometry::model::linestring<Point2D, std::vector> l;
        for (auto it = boost::geometry::ever_circling_iterator(hull.outer().begin()
                    , hull.outer().end(), false);; ++it) {
            if (std::abs(max_x_coord - (*it).x) < EPSILON
                    && ((*it).y - min_y_coord) < BOTTOM_ERR) { 
               std::copy_n(it, hull.outer().size() - 1, std::back_inserter(l));
               break;
            }
        } 

        std::string file = std::string("init_examples/contour_") + std::to_string(step);
        std::ofstream print_res(file);
        if (!print_res.is_open()) {
            std::cerr << "problem with opening : " << file << std::endl;
        }
        std::ranges::copy_if(l, std::ostream_iterator<Point2D>(print_res, "\n")
                , [min_y_coord = min_y_coord, min_x_coord =  min_x_coord
                , max_x_coord =  max_x_coord](Point2D& p) { 
                return ((p.y - min_y_coord) >= BOTTOM_ERR) || ((p.x - min_x_coord) < EPSILON) ||
                ((max_x_coord - p.x) < EPSILON); });
    }

}




