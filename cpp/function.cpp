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

void grind(std::vector<Point>& v, Point point_1, Point point_2) {
    double dist = boost::geometry::distance(point_1, point_2);
    double u_dist = dist / (POINTS_BETWEEN + 1);
    Point u_vec = (point_2 - point_1) / dist;
    for(int i : std::views::iota(1, POINTS_BETWEEN + 1)) {
       v.push_back(point_1 + i * u_vec * u_dist);
    }
}


void read_points_in_vector(std::vector<Point>& v, std::ifstream& in) {
    std::istream_iterator<Point> ii(in);
    Point point_1 = *ii;
    v.push_back(point_1);
    ++ii;

    for(Point point_2; ii != std::istream_iterator<Point>{}; ++ii) {
        point_2 = *ii;
        grind(v, point_1, point_2);
        v.push_back(point_2);
        point_1 = point_2;
    }
}

void spoil_and_get_protrusions(std::vector<Point>& v
        , std::vector<boost::geometry::model::polygon<Point, false, true, std::vector>>& vec_of_pol) {
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis1(-boost::math::constants::pi<double>()/16
            , boost::math::constants::pi<double>()/16);
    std::uniform_real_distribution<> dis2(0.03, 0.15);

    boost::geometry::model::polygon<Point, false, true, std::vector> poly_of_spoiled_points;

    int flag = -1, step = 5;
    size_t j = 0, end = v.size() - 1;
    srand(time(NULL));

    for(size_t i = 1; i < end; ++i) {
        if ((i % 40) == 1 && (i + step < end)) {
            flag = (rand()%2) ? -1 : 1;
            size_t index = i - 1;
            poly_of_spoiled_points.outer().push_back(v[i - 1]);
            for(j = i; j < (i + step); ++j) {
                Point direction = v[j] - v[j + 1];
                direction = Point(direction.y, -direction.x) * flag;
                direction = direction / boost::geometry::distance(v[j], v[j + 1]);
                direction = rotate(direction, dis1(gen));
                Point vec = direction * dis2(gen);
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


std::pair<double, double> write_intersection(std::vector<boost::geometry::model::polygon<Point>>& dif
        , std::ofstream& out_intersection) {
    double dif_area = 0, dif_perimetr = 0;
    for(auto& pol : dif) {
        std::ranges::copy(pol.outer(), std::ostream_iterator<Point>(out_intersection, "\n"));
        dif_area += boost::geometry::area(pol);
        dif_perimetr += boost::geometry::perimeter(pol);
        out_intersection << "\n\n";
    }
    return std::make_pair(dif_area, dif_perimetr);
}

