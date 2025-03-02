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
typedef Polyhedron::Face_handle                             Face_handle;

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
    void add_face(const std::map<Face_handle, norm_face>& all_faces) {
        for(const auto& val: all_faces) {
            faces.push_back(val.second);
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
        //std::cerr << el.unit_n.z() << std::endl;
    }

    size_t sz = norms.size(), start = 0, end = 0;

    struct gap_helper {
        size_t start = 0;
        double gap = 0;
    };

    std::vector<gap_helper> gaps;
    gaps.reserve(sz - 1);

    for(size_t i = 0; i < sz - 1; ++i) {
        double gap = norms[i].unit_n.z() - norms[i + 1].unit_n.z();
        gap_helper g{i, gap};
        gaps.push_back(g);
    }


    std::sort(gaps.begin(), gaps.end(), [](const gap_helper& l, const gap_helper& r) {
            return l.gap > r.gap; });

    if (gaps.size() < 2) {
        std::cerr << "Not enough elements in vector of gaps" << std::endl;
        return -1;
    }

    start = std::min(gaps[0].start, gaps[1].start);
    end = std::max(gaps[0].start, gaps[1].start);
    ++start;


    for(size_t i = start; i <= end; ++i) {
        //std::cout << norms[i].unit_n.z() << std::endl;
    }

    //std::vector<norm_face> up_rundist;
    //std::vector<norm_face> low_rundist;

    std::map<Face_handle, norm_face> up_rundist;
    std::map<Face_handle, norm_face> low_rundist;


    double z_start = norms[start].unit_n.z();
    double z_end = norms[end].unit_n.z();
    for(size_t i = start; i <= end; ++i) {
        auto begin = norms[i].it->facet_begin();
        for(size_t j = 0; j < norms[i].it->size(); ++j, ++begin) {
            auto op_face = begin->opposite()->facet();
            auto h = op_face->halfedge();
            Vector_3 unit_v = CGAL::unit_normal(h->vertex()->point(), h->next()->vertex()->point()
                    , h->next()->next()->vertex()->point());
            norm_face nf{unit_v, op_face};
            if (unit_v.z() > z_start) {
                up_rundist.insert(std::pair(op_face, nf));
            } 
            if (unit_v.z() < z_end) {
                low_rundist.insert(std::pair(op_face, nf));
            }
        }
    }
    std::cout << up_rundist.size() << std::endl;
    std::cout << low_rundist.size() << std::endl;


    Polyhedron rundist;
    Build_triangle<HalfedgeDS> triangle;
    triangle.add_face(norms, start, end);
    rundist.delegate(triangle);

    Polyhedron rundist_up;
    Build_triangle<HalfedgeDS> triangle1;
    triangle1.add_face(up_rundist);
    rundist_up.delegate(triangle1);

    Polyhedron rundist_low;
    Build_triangle<HalfedgeDS> triangle2;
    triangle2.add_face(low_rundist);
    rundist_low.delegate(triangle2);

    std::ofstream of("rundist.ply");
    of << std::fixed << std::setprecision(12);
    CGAL::IO::write_PLY(of, rundist);

    std::ofstream of1("rundist_up.ply");
    of1 << std::fixed << std::setprecision(12);
    CGAL::IO::write_PLY(of1, rundist_up);

    std::ofstream of2("rundist_low.ply");
    of2 << std::fixed << std::setprecision(12);
    CGAL::IO::write_PLY(of2, rundist_low);



    return 0;
}




    
        

