#include "../include/segment2d.h"
#include "../include/geom_functions.h"

bool Segment2D::check_line(const Segment2D& other) {
    Point2D v1 = p1 - p2;
    Point2D v2 = other.p1 - other.p2;
    if (std::abs(v1.x * v2.y - v2.x * v1.y) < EPSILON) {
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

std::ostream& operator << (std::ostream& os, const Segment2D& seg) {
    os << seg.p1 << '\n' << seg.p2;
    return os;
}


