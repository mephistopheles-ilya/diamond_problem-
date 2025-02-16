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
    std::string file_ply;
    std::string file_obj;
    boost::program_options::options_description desc("All options");
    desc.add_options()
        ("file_obj", boost::program_options::value<std::string>(&file_obj)
         -> default_value("poly.obj"), "output file to store model in obj format")
        ("file_ply", boost::program_options::value<std::string>(&file_ply)
         -> default_value("poly.ply"), "initial file with model in ply format")
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
    CGAL::IO::read_PLY(file_ply, poly);

    std::ofstream of(file_obj);
    of << std::fixed << std::setprecision(12);
    CGAL::IO::write_OBJ(of, poly);

    return 0;
}
