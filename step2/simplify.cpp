#include <CGAL/Polyline_simplification_2/Stop_above_cost_threshold.h>
#include <CGAL/Polyline_simplification_2/Stop_below_count_threshold.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/IO/STL.h>
 
 
namespace PS = CGAL::Polyline_simplification_2;
 
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                   Polygon_2;
typedef K::Point_2                           Point_2;
typedef CGAL::Polygon_with_holes_2<K>        Polygon_with_holes_2;
//typedef PS::Stop_above_cost_threshold        Stop;
//typedef PS::Stop_below_count_threshold       Stop;
typedef PS::Stop_below_count_ratio_threshold Stop;
typedef PS::Squared_distance_cost            Cost;
//typedef PS::Scaled_squared_distance_cost     Cost;
//typedef PS::Hybrid_squared_distance_cost<int>    Cost;

int main(int argc, char* argv[])
{
  std::ifstream ifs(argv[1]);
  std::istream_iterator<Point_2> it(ifs);
  Polygon_2 polygon(it, std::istream_iterator<Point_2>{});
  Polygon_with_holes_2 hpolygon(polygon);
  std::boolalpha(std::cout);
  Cost cost;
  //double a;
  //sscanf(argv[2], "%lf", &a);
  //int a = hpolygon.outer_boundary().size();
  hpolygon = PS::simplify(hpolygon, cost, Stop(0.05));
  //std::cout << hpolygon.has_holes() << std::endl;
  std::copy(hpolygon.outer_boundary().vertices_begin(), hpolygon.outer_boundary().vertices_end()
          , std::ostream_iterator<Point_2>(std::cout, "\n"));

  //std::cout.precision(12);
  //CGAL::IO::write_polygon_WKT(std::cout, polygon) << std::endl;
 
  return 0;
}
