#include <CGAL/Boolean_set_operations_2/symmetric_difference.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <list>
#include <iomanip>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef K::Point_2                                          Point_2;
typedef CGAL::Polygon_with_holes_2<K>                  Polygon_with_holes_2;

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage <file> <file> <debug>\n";
        return -1;
    }

    std::ifstream ifs1(argv[1]);
    std::istream_iterator<Point_2> it1(ifs1);
    Polygon_2 poly1(it1, std::istream_iterator<Point_2>{});
    //poly1.push_back(*(poly1.vertices_begin()));
    std::cout << poly1.size() << std::endl;

    std::ifstream ifs2(argv[2]);
    std::istream_iterator<Point_2> it2(ifs2);
    Polygon_2 poly2(it2, std::istream_iterator<Point_2>{});
    //poly2.push_back(*(poly2.vertices_begin()));
    std::cout << poly2.size() << std::endl;

    std::list<Polygon_with_holes_2> sym_diff;
    std::cout << std::boolalpha;
    std::cout << poly1.is_valid << ' ' << poly2.is_valid() << std::endl;
    std::cout << poly1.is_simple() << ' ' << poly2.is_simple() << std::endl;
    CGAL::symmetric_difference(poly1, poly1, std::back_inserter(sym_diff));
    std::cout << sym_diff.size() << std::endl;

}


    



