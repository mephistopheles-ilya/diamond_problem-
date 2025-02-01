#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <ios>

#include <boost/program_options.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/Surface_mesh/IO/PLY.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Polyhedron_3<K>                                 Polyhedron;
typedef K::Plane_3                                            Plane;
typedef K::Point_3                                            Point_3;
typedef Polyhedron::HalfedgeDS                                HalfedgeDS;

template <class HDS>
class Build_triangle : public CGAL::Modifier_base<HDS> {
public:
    typedef typename HDS::Vertex   Vertex;
    typedef typename Vertex::Point Point;

    std::vector<Point_3> points;
    std::vector<std::vector<size_t>> faces;
    size_t num_edges_loc;


    Build_triangle() {}

    void add_face(const std::vector<Point_3>& points3d, std::ifstream& in
            , std::string line, size_t num_edges,  std::string stop = "#") {
        points = points3d;
        num_edges_loc = num_edges;
        std::vector<size_t> face;
        while (line.find(stop) == std::string::npos) {
            size_t id = 0, nv = 0;
            sscanf(line.data(), "%zu %zu", &id, &nv);
            face.clear();
            face.resize(nv);
            std::getline(in, line);
            std::istringstream iss(line);
            size_t number = 0, counter = 0;
            while(counter < nv) {
                iss >> number;
                face[counter] = number;
                ++counter;
            }
            faces.push_back(std::move(face));
            std::getline(in, line);
        }
    }

    void operator()( HDS& hds) {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        unsigned int np = points.size();
        unsigned int nf = faces.size();
        B.begin_surface(np, nf, num_edges_loc, 1);
        for(unsigned int i = 0; i < np; ++i) {
            B.add_vertex(Point(points[i].x(), points[i].y(), points[i].z()));
        }

        unsigned nv = 0;
        for(unsigned int i = 0; i < nf; ++i) {
            B.begin_facet();
            nv = faces[i].size();
            std::vector<size_t>& face = faces[i];
            for(unsigned int j = 0; j < nv; ++j) {
                B.add_vertex_to_facet(face[j]);
            }
            B.end_facet();
        }
        B.end_surface();
    }
};



int main(int argc, char* argv[]) {
    std::string file_txt;
    std::string file_obj;
    std::string file_ply;
    boost::program_options::options_description desc("All options");
    desc.add_options()
        ("file_txt", boost::program_options::value<std::string>(&file_txt)
         -> default_value("InitialModels/InitialModel6"), "file with 3d model in txt format")
        ("file_obj", boost::program_options::value<std::string>(&file_obj)
         -> default_value("poly.obj"), "file to save model in obj format")
        ("file_ply", boost::program_options::value<std::string>(&file_ply)
         -> default_value("poly.ply"), "file to save model in ply format")
        ("help", "produce help message")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    std::ifstream in(file_txt);
    if(in.is_open() == false) {
        std::cerr << "Cannot open or read file" << std::endl;
        return 2;
    }

    std::string line;
    std::string stop("#");
    std::string start0("num_vertices");
    for(;;) {
        std::getline(in, line);
        if (line.find(start0) != std::string::npos)
            break;
    }
    std::getline(in, line);
    size_t num_vertices, num_facets, num_edges;
    sscanf(line.data(), "%zu %zu %zu", &num_vertices, &num_facets, &num_edges);

    std::string start1("vertices:");
    for(;;) {
        std::getline(in, line);
        if (line.find(start1) != std::string::npos)
            break;
    }
    std::getline(in, line);
    std::getline(in, line);
    std::vector<Point_3> points3d;
    while (line.find(stop) == std::string::npos) {
        size_t number = 0;
        double x = 0, y = 0, z = 0;
        sscanf(line.data(), "%zu %lf %lf %lf", &number, &x, &y, &z);
        points3d.emplace_back(x, y, z);
        std::getline(in, line);
    }

    std::string start2("indices_of_vertices");
    for(;;) {
        std::getline(in, line);
        if (line.find(start2) != std::string::npos)
            break;
    }
    std::getline(in, line);
    Polyhedron poly;
    Build_triangle<HalfedgeDS> triangle;
    triangle.add_face(points3d, in, line, num_edges);
    poly.delegate(triangle);

    std::ofstream of1(file_obj);
    of1 << std::fixed << std::setprecision(12);
    CGAL::IO::write_OBJ(of1, poly);
    std::ofstream of2(file_ply);
    of2 << std::fixed << std::setprecision(12);
    CGAL::IO::write_PLY(of2, poly);


    return 0;
}

