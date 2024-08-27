#include "../include/geom.h"
#include "../include/point2d.h"
#include "../include/point3d.h"
#include "../include/segment2d.h"
#include "../include/edge.h"

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

void projection_convex_hull(std::vector<Point3D>& arr_points3d
        , boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull, double angle) {
        boost::geometry::model::multi_point<Point2D, std::vector> mpoint;
        make_projection(angle, arr_points3d, mpoint);
        boost::geometry::convex_hull(mpoint, hull);
}

void projection_square_complexity(std::vector<Point3D>& arr_points3d
        , std::vector<Point3D>& arr_norm3d
        , std::vector<Edge>& arr_edges3d
        , boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull, double angle) {   
    std::vector<Segment2D> arr_segments2d;
    std::list<Segment2D> zero_proj2d;
    Point3D vec_proj(std::cos(angle), std::sin(angle), 0);
    for(auto& ed : arr_edges3d) {
        Point3D p1 = arr_points3d[ed.vert1_id];
        Point3D p2 = arr_points3d[ed.vert2_id];
        Point3D norm1 = arr_norm3d[ed.facet1_id];
        Point3D norm2 = arr_norm3d[ed.facet2_id];
        double dot_prod1 = vec_proj(norm1.norm());
        double dot_prod2 = vec_proj(norm2.norm());

        auto func = [sin_a = vec_proj.y, cos_a = vec_proj.x](Point3D& p) {
            double length = std::sqrt(p.x * p.x + p.y * p.y);
            double sin_b = p.y / length;
            double cos_b = p.x / length;
            double sin_a_b = sin_a * cos_b - sin_b * cos_a;
            return Point2D(length * sin_a_b, p.z);
        }; 

        if (dot_prod1 * dot_prod2 < 0) {
            Point2D proj1 = func(p1);
            Point2D proj2 = func(p2);
            arr_segments2d.push_back(Segment2D{proj1, proj2});
        }else if (dot_prod1 * dot_prod2 <= 0) {
            Point2D proj1 = func(p1);
            Point2D proj2 = func(p2);
            zero_proj2d.push_back(Segment2D{proj1, proj2});
        }
    }

    while(zero_proj2d.size()) {
        Segment2D seg = *zero_proj2d.begin();
        zero_proj2d.pop_front();
        for(auto it = zero_proj2d.begin(); it != zero_proj2d.end();) {
            if(seg.overlap(*it)) {
                seg = seg.seg_union(*it);
                zero_proj2d.erase(it);
                it = zero_proj2d.begin();
            } else {
                ++it;
            }
        }
        arr_segments2d.push_back(seg);
    }
                
    size_t all_points = arr_segments2d.size();
    Segment2D* current_seg = &arr_segments2d[0];
    hull.outer().push_back(current_seg -> p1);
    Point2D start_point = current_seg -> p2;
    for(size_t i = 0; i < all_points; ++i) {
       for(auto& seg : arr_segments2d){
           if((boost::geometry::distance(seg.p1, start_point) < EPSILON) 
                   && (current_seg != &seg)) {
               start_point = seg.p2;
               hull.outer().push_back(seg.p1);
               current_seg = &seg;
               break;
           }else if((boost::geometry::distance(seg.p2, start_point) < EPSILON) 
                   && (current_seg != &seg)) {
               start_point = seg.p1;
               hull.outer().push_back(seg.p2);  
               current_seg = &seg;
               break;
           }
       }
    }
}



