#include "header.h"

int main(int argc, char* argv[]) {
    int count = 30;
    if(argc == 2) {
        count = std::stoi(argv[1]);
    }

    for(int i = 0; i <= count; ++i) {
        std::string in_file = std::string("init_examples/contour_") + std::to_string(i);
        std::string out_file = std::string("res_examples/contour_") + std::to_string(i); 
    // FILES TO READ AND WRITE
        std::ifstream in(in_file);
        std::ofstream out(out_file);

    //LINESTRING TO KEEP INITIAL POINTS
        boost::geometry::model::linestring<Point> in_contour;
        read_points_in_contour(in_contour, in);

        std::string message;
        if(!boost::geometry::is_valid(in_contour, message)) {
            std::cout << "input contour is NOT valid " << message << std::endl;
            std::cout << "trying to correct..." << std::endl;
            boost::geometry::correct(in_contour);
            message.clear();
            if(!boost::geometry::is_valid(in_contour, message)) {
                std::cout << "input contour is still NOT valid " << message << std::endl;
            } else {
                std::cout << "input contuor is valid " << std::endl;
            }
            message.clear();
        }
        std::cout << "input linestring amount of points : " << in_contour.size() << std::endl;

    //POLYGON TO KEEP MODIFYING POINTS
        boost::geometry::model::polygon<Point, false, true, std::vector> in_polygon;
        std::ranges::for_each(in_contour, [&in_polygon](Point& p) { boost::geometry::append(in_polygon.outer(), p);});

        if(!boost::geometry::is_valid(in_polygon, message)){
            std::cout << "input polygon is NOT valid " << message << std::endl;
            std::cout << "trying to correct..." << std::endl;
            boost::geometry::correct(in_polygon);
            message.clear();
            if(!boost::geometry::is_valid(in_polygon, message)) {
                std::cout << "input polygon is still NOT valid " << message << std::endl;
            }else {
                std::cout << "input polygon is valid\n";
            }
            message.clear();
        }
        std::cout << "input polygon amount of points : " << in_polygon.outer().size() << std::endl;


    //SPOILING ANGELS OF POLYGON
        const double buffer_distance = -0.1;
        const int points_per_circle = 36;
        boost::geometry::strategy::buffer::distance_symmetric<double> distance_strategy1(buffer_distance);
        boost::geometry::strategy::buffer::join_round join_strategy1(points_per_circle);
        boost::geometry::strategy::buffer::end_round end_strategy1(points_per_circle);
        boost::geometry::strategy::buffer::point_circle circle_strategy1(points_per_circle);
        boost::geometry::strategy::buffer::side_straight side_strategy;

        std::vector<boost::geometry::model::polygon<Point, false, true, std::vector>> tmp_pol_v;
        boost::geometry::buffer(in_polygon, tmp_pol_v,
                    distance_strategy1, side_strategy,
                    join_strategy1, end_strategy1, circle_strategy1);
        in_polygon = std::move(tmp_pol_v[0]);
        tmp_pol_v.resize(0);
        boost::geometry::strategy::buffer::distance_symmetric<double> distance_strategy2(-buffer_distance);
        boost::geometry::strategy::buffer::join_round join_strategy2(60);
        boost::geometry::strategy::buffer::end_round end_strategy2(60);
        boost::geometry::strategy::buffer::point_circle circle_strategy2(60);

        boost::geometry::buffer(in_polygon, tmp_pol_v,
                    distance_strategy2, side_strategy,
                    join_strategy2, end_strategy2, circle_strategy2);
        in_polygon = std::move(tmp_pol_v[0]);
        tmp_pol_v.resize(0);

        if(!boost::geometry::is_valid(in_polygon, message)){
            std::cout << "angle-spoiled polygon is NOT valid " << message << std::endl;
            std::cout << "trying to correct..." << std::endl;
            boost::geometry::correct(in_polygon);
            message.clear();
            if(!boost::geometry::is_valid(in_polygon, message)) {
                std::cout << "it is still NOT valid " << message << std::endl;
            } else {
                std::cout << "angle-spoiled polygon is valid\n";
            }
            message.clear();
        }
        std::cout << "angle-spoiled polygon amount of points : " << in_polygon.outer().size() << std::endl;
        
    //INCREASING AMOUNT OF POINTS
        boost::geometry::model::polygon<Point, false, true, std::vector> tmp_pol;
        boost::geometry::densify(in_polygon, tmp_pol, DISTANCE_BETWEEN_POINTS);
        if(!boost::geometry::is_valid(tmp_pol, message)) {
            std::cout << "densified polygon is NOT valid " << message << std::endl;
            std::cout << "trying to correct..." << std::endl;
            boost::geometry::correct(tmp_pol);
            message.clear();
            if(!boost::geometry::is_valid(tmp_pol, message)) {
                std::cout << "densified polygon is still NOT valid " << message << std::endl;
            } else {
                std::cout << "densified polygon is valid " << std::endl;
            }
            message.clear();
        }
        in_polygon = std::move(tmp_pol);
        boost::geometry::clear(tmp_pol);
        std::cout << "densified polygon amount of points : " << in_polygon.outer().size() << std::endl;

    //CREATING SMALL SHIFTS
        create_small_shifts(in_polygon);
        if(!boost::geometry::is_valid(in_polygon, message)){
            std::cout << "small shift polygon is NOT valid " << message << std::endl;
            std::cout << "trying to correct..." << std::endl;
            boost::geometry::correct(in_polygon);
            message.clear();
            if(!boost::geometry::is_valid(in_polygon, message)) {
                std::cout << "small shift polygon is still NOT valid " << message << std::endl;
            }else {
                std::cout << "small shift polygon is valid\n";
            }
            message.clear();
        }
        std::cout << "small shift amount of points : " << in_polygon.outer().size() << std::endl;
        
    //SPOILING SOME PARTS BY SHIFTING
        std::vector<Point> v(in_polygon.outer().begin(), in_polygon.outer().end());
        std::vector<boost::geometry::model::polygon<Point, false, true, std::vector>> vec_of_protrusions;
        std::vector<bool> mask;
        spoil_and_get_protrusions(v, vec_of_protrusions, mask);

    //TRYING TO MAKE IT LOOKS MORE ROUND
        for (auto& pol : vec_of_protrusions) {
            boost::geometry::correct(pol);
            if(boost::geometry::is_valid(pol)) {
                const double buffer_distance = DISTANCE_BETWEEN_POINTS * 1.5;
                const int points_per_circle = 100;
                boost::geometry::strategy::buffer::distance_symmetric<double> distance_strategy(buffer_distance);
                boost::geometry::strategy::buffer::join_round join_strategy(points_per_circle);
                boost::geometry::strategy::buffer::end_round end_strategy(points_per_circle);
                boost::geometry::strategy::buffer::point_circle circle_strategy(points_per_circle);
                boost::geometry::strategy::buffer::side_straight side_strategy;
                boost::geometry::buffer(pol, tmp_pol_v,
                    distance_strategy, side_strategy,
                    join_strategy, end_strategy, circle_strategy);  
                boost::geometry::clear(pol);
                //boost::geometry::simplify(tmp_pol_v[0], pol, DISTANCE_BETWEEN_POINTS / 4);
                boost::geometry::convex_hull(tmp_pol_v[0], pol); 
                //pol = std::move(tmp_pol_v[0]);
                create_small_shifts(pol);
                //std::ranges::copy(pol.outer(), std::ostream_iterator<Point>(out_intersection, "\n"));
                //out_intersection << "\n\n";
                tmp_pol_v.resize(0);
            } else {
                std::cout << "NOT valid" << std::endl;
            }
        }

    //MERGING PROTRUCSIONS TO GET ONE OUT POLYGON
        boost::geometry::model::polygon<Point, false, true, std::vector> out_polygon;
        std::ranges::for_each(in_polygon.outer()
                , [&out_polygon](Point& p) { boost::geometry::append(out_polygon.outer(), p);});
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
            std::cout << "out polygon is NOT valid " << message << std::endl;
            std::cout << "trying to correct.." << std::endl;
            boost::geometry::correct(out_polygon);
            message.clear();
            if(!boost::geometry::is_valid(out_polygon, message)) {
                std::cout << "out polygon is still NOT valid " << message << std::endl;
            } else {
                std::cout << "out polygon is valid " << std::endl;
            }
            message.clear();
        }
        std::cout << "out polygon amount of points : " << out_polygon.outer().size() << std::endl;

        boost::geometry::densify(out_polygon, tmp_pol, DISTANCE_BETWEEN_POINTS);
        if(!boost::geometry::is_valid(tmp_pol, message)) {
            std::cout << "densified polygon is NOT valid " << message << std::endl;
            std::cout << "trying to correct " << std::endl;
            boost::geometry::correct(tmp_pol);
            message.clear();
            if(!boost::geometry::is_valid(tmp_pol, message)) {
                std::cout << "densified polygon is still NOT valid " << message << std::endl;
            } else {
                std::cout << "densified polygon is valid " << std::endl;
            }
            message.clear();
        }
        out_polygon = std::move(tmp_pol);
        boost::geometry::clear(tmp_pol);
        std::cout << "densifyed out polygon amount of points : " << out_polygon.outer().size() << std::endl;

    //SPOILING ANGELS OF POLYGON
        boost::geometry::buffer(out_polygon, tmp_pol_v,
                    distance_strategy1, side_strategy,
                    join_strategy1, end_strategy1, circle_strategy1);
        in_polygon = std::move(tmp_pol_v[0]);
        tmp_pol_v.resize(0);

        boost::geometry::buffer(out_polygon, tmp_pol_v,
                    distance_strategy2, side_strategy,
                    join_strategy2, end_strategy2, circle_strategy2);
        in_polygon = std::move(tmp_pol_v[0]);
        tmp_pol_v.resize(0);

        if(!boost::geometry::is_valid(out_polygon, message)){
            std::cout << "angle-spoiled polygon is NOT valid " << message << std::endl;
            std::cout << "trying to correct..." << std::endl;
            boost::geometry::correct(in_polygon);
            message.clear();
            if(!boost::geometry::is_valid(in_polygon, message)) {
                std::cout << "it is still NOT valid " << message << std::endl;
            } else {
                std::cout << "angle-spoiled polygon is valid\n";
            }
            message.clear();
        }
        std::cout << "angle-spoiled polygon amount of points : " << out_polygon.outer().size() << std::endl;

        
        std::ranges::copy(out_polygon.outer(), std::ostream_iterator<Point>(out, "\n"));
    }        
}





