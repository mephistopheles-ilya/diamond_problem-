#include <iostream>
#include <fstream>
#include <ios>

#include <boost/program_options.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/Surface_mesh/IO/PLY.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Polyhedron_3<K>                                 Polyhedron;
typedef K::Plane_3                                            Plane;
typedef K::Point_3                                            Point_3;


int main(int argc, char* argv[]) {
    std::string file_obj;
    std::string file_txt;
    std::string file_ply;
    boost::program_options::options_description desc("All options");
    desc.add_options()
        ("file_obj", boost::program_options::value<std::string>(&file_obj)
         -> default_value("poly.obj"), "file with 3d model in obj format")
        ("file_txt", boost::program_options::value<std::string>(&file_txt)
         -> default_value("poly.txt"), "output file in txt format")
        ("file_ply", boost::program_options::value<std::string>(&file_ply)
         -> default_value("poly.ply"), "output file in ply format to plot in meshlab")
        ("help", "produce help message")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    Polyhedron poly;
    CGAL::IO::read_OBJ(file_obj, poly);

    std::ofstream of(file_txt);
    of << std::fixed << std::setprecision(12);
    of << "# POLYHEDRON:" << std::endl;
    of << "# num_vertices   num_facets   num_edges" << std::endl;
    unsigned int num_vertices = poly.size_of_vertices();
    unsigned int num_facets = poly.size_of_facets();
    unsigned int num_edges = poly.size_of_halfedges() / 2;
    of << ' ' << num_vertices << ' ' << num_facets << ' ' << num_edges << std::endl;

    of << "# vertices:" << std::endl;
    of << "#   id   x y z" << std::endl;

    unsigned int i = 0;
    for(auto it = poly.vertices_begin(); it != poly.vertices_end(); ++it, ++i) {
        of << "  " << i << ' ' << (*it).point().x() << ' ' << (*it).point().y() << ' ' << (*it).point().z() << std::endl;
    }
    of << "# facets:" << std::endl;
    of << "#   id   num_of_sides   plane_coeff  (ax + by + cz + d = 0)" << std::endl;
    of << "#   indices_of_vertices" << std::endl;

    i = 0;
    for(auto fit = poly.facets_begin(); fit != poly.facets_end(); ++fit, ++i) {
        of << "   ";
        auto h = (*fit).halfedge();
        Plane plane  = Plane(h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
        double norm = std::sqrt(plane.a() * plane.a() + plane.b() * plane.b() + plane.c() * plane.c());
        auto pit = (*fit).facet_begin();
        unsigned int nv = (*fit).size();
        of << ' ' <<  i << "  " << nv << "  " << plane.a()/norm << ' ' << plane.b()/norm << ' ' <<
            plane.c()/norm << ' ' << plane.d()/norm << std::endl;
        of << "     ";
        for(unsigned int j = 0; j < nv; ++j) {
            unsigned int number = 0;
            for(auto it = poly.vertices_begin(); it != pit->vertex(); ++it, ++number) {
            }
            ++pit;
            of << number << ' ';
        }
        of << std::endl;
    }

    of << "# edges:" << std::endl;
    of << "#    vert1_id  vert2_id  facet1_id  facet2_id" << std::endl;

    for(auto edge = poly.edges_begin(); edge != poly.edges_end(); ++edge) {
        auto face1 = edge->face();
        auto face2 = edge->opposite()->face();
        auto ver1 = edge->vertex();
        auto ver2 = edge->opposite()->vertex();
        unsigned int number = 0;
        for(auto it = poly.vertices_begin(); it != ver1; ++it, ++number) {
        }
        of << number << ' ';
        number = 0;
        for(auto it = poly.vertices_begin(); it != ver2; ++it, ++number) {
        }
        of << number << ' ';
        number = 0;
        for(auto it = poly.facets_begin(); it != face1; ++it, ++number) {
        }
        of << number << ' ';
        number = 0;
        for(auto it = poly.facets_begin(); it != face2; ++it, ++number) {
        }
        of << number << std::endl;
    }
    std::ofstream of1(file_ply);
    of1 << std::fixed << std::setprecision(12);
    CGAL::IO::write_PLY(of1, poly);
    return 0;
}
