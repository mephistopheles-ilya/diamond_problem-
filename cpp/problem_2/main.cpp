#include "points.h"
#include "geom.h"

#include <vector>
#include <ranges>
#include <fstream>
#include <iomanip>
#include <numbers>
#include <random>

int main(void) {
// CREATING STREAMS ON FILES
    std::ifstream in("examples/points3d.txt");
    std::ofstream out_before("examples/points2d_before.txt");
    std::ofstream out_after("examples/points2d_after.txt");

    if(in.is_open() && out_before.is_open() && out_after.is_open()) {
        std::cout << "Files opened\n";
    } else {
        std::cout << "Wrong files\n";
        throw std::runtime_error("Cannot open or read file\n");
    }

    in >> std::fixed >> std::setprecision(PRECISION);
    out_before << std::fixed << std::setprecision(PRECISION);
    out_after << std::fixed << std::setprecision(PRECISION);


//READINF 3D POINTS IN VECTOR
    std::vector arr_points3d(std::istream_iterator<Point3D>(in), std::istream_iterator<Point3D>{});
    std::cout << "amout of points : " <<  arr_points3d.size() << '\n';


//MAKING PRIJECTION
    boost::geometry::model::multi_point<Point2D, std::vector> mpoint;
    //double radians = std::numbers::pi_v<double>/4;
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0, std::numbers::pi_v<double>);
    double radians = dis(gen);
    std::cout << "angle in radians = " << radians << std::endl;

    make_projection(radians, arr_points3d, mpoint);


//CREATING CONVEX HULL
    boost::geometry::model::polygon<Point2D, false, true, std::vector> hull;
    boost::geometry::convex_hull(mpoint, hull);


//CHECKING CONVEX HULL
    std::string message;
    if(!boost::geometry::is_valid(hull)) {
        std::cout << "convex hull is Not valid " << message;
    }else {
        std::cout << "convex hull is Valid, amount of points : " << hull.outer().size()
            << ' ' << hull.inners().size() << '\n';
    }


//FILLING POINSTS2D BEFORE
    std::ranges::copy(hull.outer(), std::ostream_iterator<Point2D>(out_before, "\n"));


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
   std::ranges::copy_if(l, std::ostream_iterator<Point2D>(out_after, "\n")
            , [min_y_coord, min_x_coord, max_x_coord](Point2D& p) { 
            return ((p.y - min_y_coord) >= BOTTOM_ERR) || ((p.x - min_x_coord) < EPSILON) ||
            ((max_x_coord - p.x) < EPSILON); });

}




