// g++-11 -std=c++23 -O2 -Wall -Wextra -Wpedantic -DBOOST_ALLOW_DEPRECATED_HEADERS main.cpp function.cpp

#include "header.h"


int main(void) {
    std::ifstream in("examples/in_decagon.txt");
    std::ofstream out("examples/out_decagon.txt");
    std::ofstream out_intersection("examples/decagon_intersection.txt");

    std::vector<Point> v;
    read_points_in_vector(v, in);
    std::cout << "Amount of points : " << v.size() << '\n';

    boost::geometry::model::polygon<Point, false, true, std::vector> polygon1;
    boost::geometry::model::polygon<Point, false, true, std::vector> polygon2;

    std::ranges::for_each(v, [&polygon1](Point& p) { boost::geometry::append(polygon1.outer(), p);});
    std::ranges::for_each(v, [&polygon2](Point& p) { boost::geometry::append(polygon2.outer(), p);});


    std::string message;
    if(!boost::geometry::is_valid(polygon1, message)){
        std::cout << "input polygon is NOT valid " << message << std::endl;
        boost::geometry::correct(polygon1);
        message.clear();
        if(!boost::geometry::is_valid(polygon1, message))
            std::cout << "it is still incorrect " << message << std::endl;
        else
            std::cout << "input polygon is correct\n";
    } else {
        std::cout << "input polygon is correct\n";
    }

    std::vector<boost::geometry::model::polygon<Point, false, true, std::vector>> vec_of_protrusions;
    spoil_and_get_protrusions(v, vec_of_protrusions);

    for (auto& pol : vec_of_protrusions) {
        boost::geometry::correct(pol);
        if(boost::geometry::is_valid(pol)) {
			boost::geometry::model::polygon<Point, false, true, std::vector> hull;
            boost::geometry::convex_hull(pol, hull);
            std::ranges::copy(hull.outer(), std::ostream_iterator<Point>(out_intersection, "\n"));
            out_intersection << "\n\n";
        } else {
            std::cout << "NOT valid" << std::endl;
        }
    }
        
    std::ranges::copy(v, std::ostream_iterator<Point>(out, "\n"));

    double area = 0, perimetr = 0;
    for (auto& pol : vec_of_protrusions) {
        if(boost::geometry::is_valid(pol)) {
            area += boost::geometry::area(pol);
            perimetr += boost::geometry::perimeter(pol);
        }
    }

    std::cout << "arae : " << area << std::endl;
    std::cout << "perimetr : " << perimetr << std::endl;
    std::cout << "ratio : " << area / perimetr << std::endl;
}









