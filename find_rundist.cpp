#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <map>

#include <boost/program_options.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K>                               Polyhedron;
typedef K::Point_3                                          Point_3;
typedef K::Plane_3                                          Plane;
typedef K::Vector_3                                         Vector_3;
typedef Polyhedron::Vertex_handle                           Vertex_handle;
typedef Polyhedron::HalfedgeDS                              HalfedgeDS;

struct norm_face {
    Vector_3 unit_n;
    typename Polyhedron::Face_iterator it;
};


template <class HDS>
class Build_triangle : public CGAL::Modifier_base<HDS> {
public:
    typedef typename HDS::Vertex   Vertex;
    typedef typename Vertex::Point Point;

    std::vector<norm_face> faces;


    Build_triangle() {}

    void add_face(const std::vector<norm_face>& all_faces, size_t begin, size_t end) {
        for(size_t i = begin; i <= end; ++i) {
            faces.push_back(all_faces[i]);
        }
    }

    void operator()( HDS& hds) {
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        std::map<Vertex_handle, size_t> vertex_map;
        size_t vertex_index = 0;
        for(auto& face: faces) {
            auto begin = face.it->facet_begin();
            size_t num = face.it->size();
            for(size_t i = 0; i < num; ++i, ++begin) {
                if (vertex_map.find(begin->vertex()) == vertex_map.end()) {
                    vertex_map[begin->vertex()] = vertex_index;
                    ++vertex_index;
                }
            }
        }
        B.begin_surface(vertex_map.size(), faces.size(), vertex_map.size(), 1);
        std::vector<std::pair<Vertex_handle, size_t>> v;
        std::copy(vertex_map.begin(), vertex_map.end(), std::back_inserter(v));
        std::sort(v.begin(), v.end(), [](const std::pair<Vertex_handle, size_t>& a
                    , const std::pair<Vertex_handle, size_t>& b) { return a.second < b.second;});
        for(const auto& el: v) {
            B.add_vertex(el.first->point());
        }


        for(auto& face: faces) {
            B.begin_facet();
            auto begin = face.it->facet_begin();
            size_t num = face.it->size();
            for(size_t i = 0; i < num; ++i, ++begin) {
                B.add_vertex_to_facet(vertex_map[begin->vertex()]);
            }
            B.end_facet();
        }
        B.end_surface();
    }
};


int main(int argc, char* argv[]) {
    std::string type_of_file;
    std::string file;
    boost::program_options::options_description desc("All options");
    desc.add_options()
        ("type_of_file", boost::program_options::value<std::string>(&type_of_file)
         -> default_value("ply"), "ply or obj")
        ("file", boost::program_options::value<std::string>(&file)
         -> default_value("poly.ply"), "first polyhedron")
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
    if (type_of_file == "ply") {
        CGAL::IO::read_PLY(file, poly);
    } else {
        CGAL::IO::read_OBJ(file, poly);
    }

    std::vector<norm_face> norms;
    for(auto fit = poly.facets_begin(); fit != poly.facets_end(); ++fit) {
        auto h = (*fit).halfedge();
        Vector_3 unit_v = CGAL::unit_normal(h->vertex()->point(), h->next()->vertex()->point()
                , h->next()->next()->vertex()->point());
        norms.emplace_back(unit_v, fit);
    }

    std::sort(norms.begin(), norms.end(), [](const norm_face& n1, const norm_face& n2)
            { return n1.unit_n.z() > n2.unit_n.z();});
    for(auto& el: norms) {
        std::cerr << el.unit_n.z() << std::endl;
    }
    double max_pl_dis = 0, pl_dis = 0;
    double max_mi_dis = 0, mi_dis = 0;
    size_t sz = norms.size(), start = 0, end = 0;
    for(size_t i = 0; i < sz - 1; ++i) {
        if (norms[i + 1].unit_n.z() > 0) {
            pl_dis = norms[i].unit_n.z() - norms[i + 1].unit_n.z();
            if (pl_dis > max_pl_dis) {
                max_pl_dis = pl_dis;
                start = i + 1;
            }
        }
        if (norms[i].unit_n.z() < 0) {
            mi_dis = norms[i].unit_n.z() - norms[i + 1].unit_n.z();
            if (mi_dis > max_mi_dis) {
                max_mi_dis = mi_dis;
                end = i;
            }
        }
    }
    for(size_t i = start; i <= end; ++i) {
        std::cout << norms[i].unit_n.z() << std::endl;
    }
    Polyhedron rundist;
    Build_triangle<HalfedgeDS> triangle;
    triangle.add_face(norms, start, end);
    rundist.delegate(triangle);

    std::ofstream of("rundist.ply");
    of << std::fixed << std::setprecision(12);
    CGAL::IO::write_PLY(of, rundist);


    return 0;
}




    
        

