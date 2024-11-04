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


#include "include/point2d.hpp"
#include "include/geom_functions.hpp"

using Polygon_2 = boost::geometry::model::polygon<Point2D, false, true, std::vector>;
 
int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "First and second arguments are file\n";
        return 0;
    }
    std::ifstream in1(argv[1]);
    std::istream_iterator<Point2D> it1(in1);
    Polygon_2 poly1;
    std::copy(it1, std::istream_iterator<Point2D>{}, std::back_inserter(poly1.outer()));
    poly1.outer().push_back(*(poly1.outer().begin()));

    std::ifstream in2(argv[2]);
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
    boost::geometry::sym_difference(poly1, poly2, polys);
    std::cout << "amount of polygons : " << polys.size() << std::endl;
    double area = 0;
    for(auto& pol: polys) {
        //for(auto it = pol.outer().begin(); it != pol.outer().end(); ++it) {
        //    std::cout << *it << std::endl;
        //}
        area += boost::geometry::area(pol);
        //std::cout << std::endl << std::endl;
    }
    double per = boost::geometry::perimeter(poly1) + boost::geometry::perimeter(poly2);
    per /= 2;
    std::cout << area / per << std::endl;
}

