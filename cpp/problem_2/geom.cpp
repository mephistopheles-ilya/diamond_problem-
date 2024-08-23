#include "geom.h"
#include "points.h"

void make_projection(double angle, std::vector<Point3D>& arr_points3d
        , boost::geometry::model::multi_point<Point2D, std::vector>& mpoint) {
    double sin_a = std::sin(angle);
    double cos_a = std::cos(angle);
    auto func = [sin_a, cos_a](Point3D& p) {
        double length = std::sqrt(p.x * p.x + p.y * p.y);
        double sin_b = p.y / length;
        double cos_b = p.x / length;
        double sin_a_b = sin_a * cos_b - sin_b * cos_a;
        return Point2D(length * sin_a_b, p.z);
    }; 
    std::ranges::transform(arr_points3d, std::back_inserter(mpoint), func);
}

Point2D find_minimal_y(boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull) {
    Point2D p_min_y = hull.outer()[0];
    for (auto& p : hull.outer()) {
        if (p.y < p_min_y.y) 
            p_min_y = p;
    }
    return p_min_y;
}

std::pair<double, double> find_borders_for_x(Point2D lowest
        , boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull) {
    double min_y_coord = lowest.y;
    double min_x_coord = 100;
    double max_x_coord = -100;
    for (auto& p : hull.outer()) {
        if ((p.y - min_y_coord) < BOTTOM_ERR) {
            if (p.x < min_x_coord)
                min_x_coord = p.x;
            if (p.x > max_x_coord)
                max_x_coord = p.x;
        }
    }
    return std::make_pair(min_x_coord, max_x_coord);
}


