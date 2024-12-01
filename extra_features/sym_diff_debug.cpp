#include <boost/geometry/algorithms/detail/is_valid/interface.hpp>
#include <boost/geometry/algorithms/perimeter.hpp>
#include <boost/geometry/algorithms/sym_difference.hpp>
#include <iostream>
#include <iterator>
#include <boost/geometry.hpp>
#include <ostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <boost/geometry/geometries/polygon.hpp> 
#include <ranges>


#include "include/point2d.hpp"
#include "include/geom_functions.hpp"

using Polygon_2 = boost::geometry::model::polygon<Point2D, false, true, std::vector>;
 
int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage <file> <file> <debug>\n";
        return -1;
    }

    std::string contour1 = std::string(argv[1]);
    std::ifstream in1(contour1);
    std::istream_iterator<Point2D> it1(in1);
    Polygon_2 poly1;
    std::copy(it1, std::istream_iterator<Point2D>{}, std::back_inserter(poly1.outer()));
    poly1.outer().push_back(*(poly1.outer().begin()));

    std::string contour2 = std::string(argv[2]);
    std::ifstream in2(contour2);
    std::istream_iterator<Point2D> it2(in2);
    Polygon_2 poly2;
    std::copy(it2, std::istream_iterator<Point2D>{}, std::back_inserter(poly2.outer()));
    poly2.outer().push_back(*(poly2.outer().begin()));

    if(boost::geometry::is_valid(poly1) != true) {
        std::cout << "poly1 is not valid" << std::endl;
    }
    if(boost::geometry::is_valid(poly2) != true) {
        std::cout << "poly2 is not valid" << std::endl;
    }

    std::vector<Polygon_2> polys;
    boost::geometry::sym_difference(poly2, poly1, polys);
    std::ofstream out(argv[3]);
    std::ostream_iterator<Point2D> oit(out, "\n");
    double area = 0;
    std::cout << "Amount of polygons : " << polys.size() << std::endl;
    for(auto& pol: polys) {
        area += boost::geometry::area(pol);
        std::ranges::copy(pol.outer(), oit);
        out << *(pol.outer().begin()) << std::endl << std::endl;
    }
    double per = boost::geometry::perimeter(poly1) + boost::geometry::perimeter(poly2);
    per /= 2;
    std::cout << "Difference : " << area / per << std::endl;
}

