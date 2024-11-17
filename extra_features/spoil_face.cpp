#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Polyhedron_3_fwd.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/IO/STL.h>
#include <CGAL/IO/PLY.h>

#include <cmath>
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
    srand(time(NULL));
    while(line.find(stop) == std::string::npos) {
        double a, b, c, d;
        int unused1, unused2;
        sscanf(line.c_str(), "%d %d %lf %lf %lf %lf", &unused1, &unused2, &a, &b, &c, &d);
        if (c > 0) {
            typename K::Plane_3 plane(a, b, c, d);
            planes.push_back(plane);
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
    std::list<Plane> planes;
    std::ifstream is(argv[1]);
    read_planes(is, planes);
 
    Polyhedron_3 chull;

    CGAL::halfspace_intersection_with_constructions_3(planes.begin(), planes.end(), chull);

#if 0
    std::ofstream of1("tmp.stl");
    //std::ofstream of2("tmp.ply");
    //CGAL::IO::write_STL(of1, chull);
    //CGAL::IO::write_PLY(of2, chull);

    of1 << "solid\n";
    for(auto face = chull.facets_begin(); face != chull.facets_end(); ++face) {
        auto it_point = face->facet_begin();
        Point p1 = it_point->vertex()->point();
        ++it_point;
        Point p2 = it_point->vertex()->point();
        ++it_point;
        Point p3 = it_point->vertex()->point();
        typename K::Plane_3 plane(p1, p2, p3);
        double length = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
        of1 << "facet normal " << plane.a()/length << ' ' << plane.b()/length << ' ' 
            << plane.c()/length << "\nouter loop\n";
        of1 << "vertex " << p1 << "\n";
        of1 << "vertex " << p2 << "\n";
        of1 << "vertex " << p3 << "\n";

        for(int i = 3; i < face->size(); ++i) {
            ++it_point;
            of1 << "vertex " << it_point->vertex()->point() << "\n";
        }
        of1 << "endloop\nendfacet\n";
#endif

#if 0
        for(int i = 3; i < face->size(); ++i) {
            of1 << "facet normal " << plane.a()/length << ' ' << plane.b()/length << ' ' 
                    << plane.c()/length << "\nouter loop\n";
            p2 = p3;
            p3 = (++it_point)->vertex()->point();
            of1 << "vertex " << p1 << "\n";
            of1 << "vertex " << p2 << "\n";
            of1 << "vertex " << p3 << "\n";
            of1 << "endloop\nendfacet\n";
        }
    }
    of1 << "endsolid"<<std::endl;
#endif

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

#if 0
    int counter = 0;
    for(auto it = chull.edges_begin(); it != chull.edges_end(); ++it) {
        std::cout << it->vertex()->point() << std::endl;
        std::cout << it->opposite()->vertex()->point() << std::endl;
        ++counter;
        if (counter == 10) break;
    }
#endif


    return 0;
}
