#include "../include/segment2d.hpp"
#include "../include/geom_functions.hpp"
#include <boost/geometry/algorithms/detail/distance/interface.hpp>
#include <utility>

bool Segment2D::check_line(const Segment2D& other) {
    Point2D v1 = p1 - p2;
    Point2D v2 = other.p1 - other.p2;
    double a1 = -v1.y;
    double b1 = v1.x;
    double c1 = v1.x * p1.y - p1.x * v1.y;
    if (std::abs(v1.x * v2.y - v2.x * v1.y) < EPSILON
            && std::abs(a1 * other.p1.x + b1 * other.p1.y - c1) < EPSILON) {
        return true;
    }else {
        return false;
    }
}

bool Segment2D::overlap(const Segment2D& other) {
    if (!check_line(other)) {
        return false;
    }else {
        if(p1.x <= std::max(other.p1.x, other.p2.x) 
                && p1.x >= std::min(other.p1.x, other.p2.x)) {
            return true;
        }else if(p2.x <= std::max(other.p1.x, other.p2.x) 
                && p2.x >= std::min(other.p1.x, other.p2.x)) {
            return true;
        }else if(p1.y <= std::max(other.p1.y, other.p2.y) 
                && p1.y >= std::min(other.p1.y, other.p2.y)) {
            return true;
        }else if(p2.y <= std::max(other.p1.y, other.p2.y) 
                && p2.y >= std::min(other.p1.y, other.p2.y)){
            return true;
        }else {
            return false;
        }
    }
}

Segment2D Segment2D::seg_union(const Segment2D& other) {
    Point2D pp1(std::max({p1.x, p2.x, other.p1.x, other.p2.x})
                , std::max({p1.y, p2.y, other.p1.y, other.p2.y}));
    Point2D pp2(std::min({p1.x, p2.x, other.p1.x, other.p2.x})
                , std::min({p1.y, p2.y, other.p1.y, other.p2.y}));
    return Segment2D{pp1, pp2};
}

std::pair<bool, Point2D> Segment2D::intersect(const Segment2D& other) {
    /* 
     * a1*x + b1*y = c1
     * a2*x + b2*y = c2
     */
    Point2D vec1 = p1 - p2;
    Point2D vec2 = other.p1 - other.p2;
    vec1 = vec1 / boost::geometry::distance(p1, p2);
    vec2 = vec2 / boost::geometry::distance(other.p1, other.p2);
    double a1 = -vec1.y;
    double b1 = vec1.x;
    double c1 = vec1.x * p1.y - p1.x * vec1.y;
    double a2 = -vec2.y;
    double b2 = vec2.x;
    double c2 = vec2.x * other.p1.y - other.p1.x * vec2.y;
    double det = a1 * b2 - a2 * b1;
    if(std::abs(det) >= EPSILON) {
        Point2D point = Point2D(c1 * b2 - b1 * c2, a1 * c2 - a2 * c1) / det;
        bool check = point.x <= std::max(other.p1.x, other.p2.x) 
                && point.x >= std::min(other.p1.x, other.p2.x)
                && point.y <= std::max(other.p1.y, other.p2.y) 
                && point.y >= std::min(other.p1.y, other.p2.y)
                && point.x <= std::max(p1.x, p2.x) 
                && point.x >= std::min(p1.x, p2.x)
                && point.y <= std::max(p1.y, p2.y) 
                && point.y >= std::min(p1.y, p2.y);
           
        if(check) {
            return std::make_pair(true, point);
        }else {
            return std::make_pair(false, point);
        }
    }
    return std::make_pair(false, p1);
}



std::ostream& operator << (std::ostream& os, const Segment2D& seg) {
    os << seg.p1 << '\n' << seg.p2;
    return os;
}


