#include "include/point2d.h"
#include "include/point3d.h"
#include "include/segment2d.h"
#include "include/edge.h"
#include "include/functions.h"
#include "include/geom.h"

#include <fstream>
#include <iomanip>
#include <numbers>

int main(int argc, char* argv[]) {
    int projections = 30;
    std::string file = "InitialModel";
    if (argc == 3) {
        projections = std::stoi(argv[1]);
        file = argv[2];
    }

//CREATING STREAMS ON FILES
    std::ifstream in(file);
    if(in.is_open()) {
        std::cout << "Input files are opened\n";
    } else {
        std::cout << "Wrong files\n";
        throw std::runtime_error("Cannot open or read file\n");
    }
    in >> std::fixed >> std::setprecision(PRECISION);

//READIND DATA FROM FILE
    std::vector<Point3D> arr_points3d;
    std::vector<Point3D> arr_norm3d;
    std::vector<Edge> arr_edges3d;
    auto [num_vertices, num_facest, num_edges] = read_structures_from_file(in
            , arr_points3d, arr_norm3d, arr_edges3d); 
    std::cout << "num_vertices = " << num_vertices << ' ' << arr_points3d.size() << '\n';
    std::cout << "num_facest = " << num_facest << ' ' << arr_norm3d.size() << '\n';
    std::cout << "num_edges = " << num_edges << ' ' << arr_edges3d.size() <<  std::endl;

//MAKING PRIJECTION
    for (int step = 0; step <= projections; ++step) { 
        boost::geometry::model::polygon<Point2D, false, true, std::vector> hull;
        double angle = step * std::numbers::pi_v<double>/projections;
        //projection_convex_hull(arr_points3d, hull, angle);
        projection_square_complexity(arr_points3d, arr_norm3d, arr_edges3d, hull, angle);

//CHECKING CONVEX HULL
        boost::geometry::correct(hull);
        std::string message;
        if(!boost::geometry::is_valid(hull)) {
            std::cout << "convex hull is Not valid " << message;
        }else {
            std::cout << "convex hull is Valid, amount of points : " << hull.outer().size()
                << ' ' << hull.inners().size() << '\n';
        }

//FINDING POINTS WITH MINIMAL Y COORD
        Point2D p_min_y = find_minimal_y(hull);
        std::cout << "The lowest point : x = " << p_min_y.x << " y = " << p_min_y.y << std::endl;

//FINDING POINTS FITH MAX AND MIN X CORD AND ALSO IN THE BOTTOM
        double min_y_coord = p_min_y.y;
        auto [min_x_coord, max_x_coord] = find_borders_for_x(p_min_y, hull);
        std::cout << "min y : " << min_y_coord << std::endl;
        std::cout << "max x : " << max_x_coord << std::endl;
        std::cout << "min x : " << min_x_coord << std::endl;


//FINDING BEGINING OF THE COUNTOR
        boost::geometry::model::linestring<Point2D, std::vector> l;
        for (auto it = boost::geometry::ever_circling_iterator(hull.outer().begin()
                    , hull.outer().end(), false);; ++it) {
            if (max_x_coord - (*it).x < EPSILON) { 
               std::copy_n(it, hull.outer().size() - 1, std::back_inserter(l));
               break;
            }
        } 

//PRINTING IN FILI RIGHT POINTS
        std::string file = std::string("init_examples/contour_") + std::to_string(step);
        std::ofstream print_res(file);
        if (print_res.is_open())
            std::cout << "file is opened\n";
        else 
            std::cout << "problem with opening\n";
        std::ranges::copy_if(l, std::ostream_iterator<Point2D>(print_res, "\n")
                , [min_y_coord, min_x_coord, max_x_coord](Point2D& p) { 
                return ((p.y - min_y_coord) >= BOTTOM_ERR) || ((p.x - min_x_coord) < EPSILON) ||
                ((max_x_coord - p.x) < EPSILON); });
    }
}




