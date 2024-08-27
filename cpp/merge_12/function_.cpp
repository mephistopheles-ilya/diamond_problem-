#include "header.h"

std::istream& operator >> (std::istream& is, Point& p) {
    is >> p.x >> p.y;
    return is;
}


std::ostream& operator << (std::ostream& os, const Point& p) {
    os << p.x << ' ' << p.y; 
    return os;
}

Point rotate(Point p, double angle) {
    boost::numeric::ublas::matrix<double, boost::numeric::ublas::row_major
        , std::vector<double>> rot_matrix(2, 2);
    boost::numeric::ublas::vector<double, std::vector<double>> vec(2);
    vec(0) = p.x;
    vec(1) = p.y;
    rot_matrix(0, 0) = cos(angle);
    rot_matrix(0, 1) = -sin(angle);
    rot_matrix(1, 0) = sin(angle);
    rot_matrix(1, 1) = cos(angle);
    auto res = boost::numeric::ublas::prod(rot_matrix, vec);
    return Point(res(0), res(1));
}


void read_points_in_contour(boost::geometry::model::linestring<Point>& in_contour, std::ifstream& in) {
    std::istream_iterator<Point> ii(in);
    Point point_1 = *ii;
    in_contour.push_back(point_1);
    ++ii;

    for(Point point_2; ii != std::istream_iterator<Point>{}; ++ii) {
        point_2 = *ii;
        double dist = boost::geometry::distance(point_1, point_2);
        if (dist >= EPSILON) {
        	in_contour.push_back(point_2);
        }
        point_1 = point_2;
    }
}

void create_small_shifts(boost::geometry::model::polygon<Point, false, true, std::vector>& pol) {
    auto& vec_of_points = pol.outer();
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(-DISTANCE_BETWEEN_POINTS / 20, DISTANCE_BETWEEN_POINTS / 20);
    size_t end = vec_of_points.size() - 1;
    for(size_t i = 1; i < end; ++i) {
        Point direction = vec_of_points[i] - vec_of_points[i + 1];
        direction = Point(direction.y, -direction.x);
        direction = direction / boost::geometry::distance(vec_of_points[i], vec_of_points[i + 1]);
        direction = direction * dis(gen);
        vec_of_points[i] = vec_of_points[i] + direction;
    }
}


void spoil_and_get_protrusions(std::vector<Point>& v
        , std::vector<boost::geometry::model::polygon<Point, false, true, std::vector>>& vec_of_pol
        , std::vector<bool>& mask) {
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis1(-boost::math::constants::pi<double>()/16
            , boost::math::constants::pi<double>()/16);
    std::uniform_real_distribution<> dis_small(DISTANCE_BETWEEN_POINTS/2 , DISTANCE_BETWEEN_POINTS);
    std::uniform_real_distribution<> dis_bigger(DISTANCE_BETWEEN_POINTS, DISTANCE_BETWEEN_POINTS * 2);

    boost::geometry::model::polygon<Point, false, true, std::vector> poly_of_spoiled_points;

    int flag = -1, step = 8;
    size_t j = 0, end = v.size() - 20;
    srand(time(NULL));

    for(size_t i = 1; i < end; ++i) {
        if ((i % 197) == 30 && (i + step < end) && (rand()%3 == 1)) {
            flag = (rand()%2) ? -1 : 1;
            if (flag == 1) mask.push_back(1);
            else mask.push_back(0);
            size_t index = i - 1;
            poly_of_spoiled_points.outer().push_back(v[i - 1]);
            for(j = i; j < (i + step); ++j) {
                Point direction = v[j] - v[j + 1];
                direction = Point(direction.y, -direction.x) * flag;
                direction = direction / boost::geometry::distance(v[j], v[j + 1]);
                //direction = rotate(direction, dis1(gen));
                Point vec = direction * (((j - i) < 3 || (j - i) > 5 ) ? dis_small(gen) : dis_bigger(gen));
                v[j] = v[j] + vec;
                poly_of_spoiled_points.outer().push_back(v[j]);
            }
            poly_of_spoiled_points.outer().push_back(v[j + 1]);
            poly_of_spoiled_points.outer().push_back(v[index]);
            vec_of_pol.push_back(poly_of_spoiled_points);
            boost::geometry::clear(poly_of_spoiled_points);
            i = j - 1;
       }
    }
}


std::pair<double, double> write_intersection(std::vector<boost::geometry::model::polygon
        <Point, false, true, std::vector>>& dif, std::ofstream& out_intersection) {
    double dif_area = 0, dif_perimetr = 0;
    for(auto& pol : dif) {
        std::ranges::copy(pol.outer(), std::ostream_iterator<Point>(out_intersection, "\n"));
        dif_area += boost::geometry::area(pol);
        dif_perimetr += boost::geometry::perimeter(pol);
        out_intersection << "\n\n";
    }
    return std::make_pair(dif_area, dif_perimetr);
}

