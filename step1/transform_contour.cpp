#include "include/point2d.h"
#include "include/geom_functions.h"
#include "include/rw_functions.h"

#include <algorithm>
#include <boost/geometry/iterators/ever_circling_iterator.hpp>
#include <fstream>

#include <boost/program_options.hpp>
#include <iterator>

int main(int argc, char* argv[]) {
    int projections;
    bool grid;
    std::string grid_method;
    bool shift;
    bool big_shift;
    double dist_points;
    double sz_shift;
    boost::program_options::options_description desc("All options");
    desc.add_options()
        ("projections", boost::program_options::value<int>(&projections) -> default_value(30)
         , "how many projections do you want to create")
        ("grid", boost::program_options::value<bool>(&grid) 
         -> default_value(true), "do you want to increase amount of points")
        ("gr_method", boost::program_options::value<std::string>(&grid_method)
         -> default_value("uniform"), "uniform or ceils")
        ("shift", boost::program_options::value<bool>(&shift) 
         -> default_value(true), "do you want to create small shifts")
        ("big_shifts", boost::program_options::value<bool>(&big_shift) 
         -> default_value(false), "do you want to big_shifts")
        ("dist_points", boost::program_options::value<double>(&dist_points) 
         -> default_value(__DISTANCE_BETWEEN_POINTS)
         , "distnce between points after grid or size of ceils")
        ("sz_shift", boost::program_options::value<double>(&sz_shift) 
         -> default_value(0.001), "size of small shifts")
        ("help", "produce help message")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }

    for(int i = 0; i <= projections; ++i) {
        std::string in_file = std::string("init_examples/contour_") + std::to_string(i);
        std::string out_file = std::string("res_examples/contour_") + std::to_string(i); 
    // FILES TO READ AND WRITE
        std::ifstream in(in_file);
        std::ofstream out(out_file);

    //LINESTRING TO KEEP INITIAL POINTS
        boost::geometry::model::linestring<Point2D> in_contour;
        read_points_in_contour(in_contour, in);

        std::string message;
        if(!boost::geometry::is_valid(in_contour, message)) {
            std::cerr << "input contour is NOT valid " << message << std::endl;
            std::cerr << "trying to correct..." << std::endl;
            boost::geometry::correct(in_contour);
            message.clear();
            if(!boost::geometry::is_valid(in_contour, message)) {
                std::cerr << "input contour is still NOT valid " << message << std::endl;
            }
            message.clear();
        }

        boost::geometry::model::linestring<Point2D> tmp_contour;
        if (grid) {
            if(grid_method == "uniform") {
                boost::geometry::densify(in_contour, tmp_contour, dist_points);
            } else if (grid_method == "ceils") {
                uniform_grid_intersection(in_contour, tmp_contour, dist_points);
            }
            if(!boost::geometry::is_valid(tmp_contour, message)) {
                std::cout << "densified contour is NOT valid " << message << std::endl;
                std::cout << "trying to correct..." << std::endl;
                boost::geometry::correct(tmp_contour);
                message.clear();
                if(!boost::geometry::is_valid(tmp_contour, message)) {
                    std::cout << "densified contour is still NOT valid " << message << std::endl;
                } else {
                    std::cout << "densified contour is valid " << std::endl;
                }
                message.clear();
            }
            in_contour = std::move(tmp_contour);
            boost::geometry::clear(tmp_contour);
        }

        if (shift) {
            create_small_shifts(in_contour, sz_shift);
            if(!boost::geometry::is_valid(in_contour, message)){
                std::cerr << "small shift contour is NOT valid " << message << std::endl;
                std::cerr << "trying to correct..." << std::endl;
                boost::geometry::correct(in_contour);
                message.clear();
                if(!boost::geometry::is_valid(in_contour, message)) {
                    std::cerr << "small shift contour is still NOT valid " << message << std::endl;
                }
                message.clear();
            }
        }

        if (big_shift) {
            Point2D start = *in_contour.begin();
            boost::geometry::model::polygon<Point2D, false, true, std::vector> in_polygon;
            std::ranges::copy(in_contour, std::back_inserter(in_polygon.outer()));

            std::vector<Point2D> v(in_polygon.outer().begin(), in_polygon.outer().end());
            std::vector<boost::geometry::model::polygon<Point2D, false, true, std::vector>>
                vec_of_protrusions;
            std::vector<boost::geometry::model::polygon<Point2D, false, true, std::vector>>
                tmp_pol_v;
            std::vector<bool> mask;
            spoil_and_get_protrusions(v, vec_of_protrusions, mask);

            for (auto& pol : vec_of_protrusions) {
                boost::geometry::correct(pol);
                if(boost::geometry::is_valid(pol)) {
                    const double buffer_distance = __DISTANCE_BETWEEN_POINTS * 1.5;
                    const int points_per_circle = 100;
                    boost::geometry::strategy::buffer::distance_symmetric<double> 
                        distance_strategy(buffer_distance);
                    boost::geometry::strategy::buffer::join_round join_strategy(points_per_circle);
                    boost::geometry::strategy::buffer::end_round end_strategy(points_per_circle);
                    boost::geometry::strategy::buffer::point_circle 
                        circle_strategy(points_per_circle);
                    boost::geometry::strategy::buffer::side_straight side_strategy;
                    boost::geometry::buffer(pol, tmp_pol_v,
                        distance_strategy, side_strategy,
                        join_strategy, end_strategy, circle_strategy);  
                    boost::geometry::clear(pol);
                    boost::geometry::convex_hull(tmp_pol_v[0], pol); 
                    tmp_pol_v.resize(0);
                }
            }

            boost::geometry::model::polygon<Point2D, false, true, std::vector> out_polygon;
            std::ranges::for_each(in_polygon.outer()
                    , [&out_polygon](Point2D& p) {
                    boost::geometry::append(out_polygon.outer(), p);});
            int counter = 0;
            for (auto& pol : vec_of_protrusions) {
                if(boost::geometry::is_valid(pol)) {
                    if(!mask[counter]) {
                        boost::geometry::union_(out_polygon, pol, tmp_pol_v);
                        out_polygon = std::move(tmp_pol_v[0]);
                        tmp_pol_v.resize(0);
                    } else {
                        boost::geometry::difference(out_polygon, pol, tmp_pol_v);
                        out_polygon = std::move(tmp_pol_v[0]);
                        tmp_pol_v.resize(0);
                    }
                }
                ++counter;
            }
            if(!boost::geometry::is_valid(out_polygon, message)) {
                std::cerr << "out polygon is NOT valid " << message << std::endl;
                std::cerr << "trying to correct.." << std::endl;
                boost::geometry::correct(out_polygon);
                message.clear();
                if(!boost::geometry::is_valid(out_polygon, message)) {
                    std::cerr << "out polygon is still NOT valid " << message << std::endl;
                }
                message.clear();
            }
       
            for (auto it = boost::geometry::ever_circling_iterator(out_polygon.outer().begin()
                        , out_polygon.outer().end(), false);;++it) {
                    if (boost::geometry::distance(*it, start) < EPSILON) {
                        std::ranges::copy_n(it, out_polygon.outer().size() - 1
                            , std::ostream_iterator<Point2D>(out, "\n"));
                        break;
                    }
            }
        } else {
            std::ranges::copy(in_contour, std::ostream_iterator<Point2D>(out, "\n"));
        }
    }        
}





