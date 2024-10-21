#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Polyhedron_3_fwd.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/IO/STL.h>
#include <CGAL/IO/PLY.h>

#include <cstdio>
#include <list>
#include <string>
#include <istream>



typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Polyhedron_3<K>                                 Polyhedron_3;
typedef K::Plane_3                                            Plane;
typedef K::Point_3                                            Point;
typedef CGAL::Surface_mesh<Point>                             Surface_mesh;


void read_planes(std::istream& is, std::list<Plane>& planes) {
    std::string line;
    std::string key("# plane coeff");
    std::string stop("#");
    for(;;) {
        std::getline(is, line);
        if (line.find(key) != std::string::npos) 
            break;
    }
    std::getline(is, line);
    int counter = 1;
    while(line.find(stop) == std::string::npos) {
        std::getline(is, line);
        double a, b, c, d;
        sscanf(line.c_str(), "%lf %lf %lf %lf", &a, &b, &c, &d);
        if (counter == 0) {
            typename K::Plane_3 plane(a + 1, b, c, d);
            planes.push_back(plane);
        } else {
            typename K::Plane_3 plane(a, b, c, d);
            planes.push_back(plane);
        }
        std::getline(is, line);
        ++counter;
    }
}
 
int main (int argc, char* argv[]) {
    std::list<Plane> planes;
    std::ifstream is(argv[1]);
    read_planes(is, planes);
 
    Polyhedron_3 chull;

    CGAL::halfspace_intersection_with_constructions_3(planes.begin(), planes.end(), chull);

    std::ofstream of1("tmp.stl");
    std::ofstream of2("tmp.ply");
    CGAL::IO::write_STL(of1, chull);
    CGAL::IO::write_PLY(of2, chull);

#if 1
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
