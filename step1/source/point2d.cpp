#include "../include/point2d.h"
#include "../include/geom_functions.h"

Point2D operator + (const Point2D& left, const Point2D& right) noexcept {
    return Point2D(left.x + right.x, left.y + right.y);
}

Point2D operator - (const Point2D& left, const Point2D& right) noexcept {
    return Point2D(left.x - right.x, left.y - right.y);
}

Point2D operator * (const Point2D& left, double right) noexcept {
    return Point2D(left.x * right, left.y * right);
}

Point2D operator * (double left, const Point2D& right) noexcept {
    return Point2D(right.x * left, right.y * left);
}

Point2D operator / (const Point2D& left, double right) noexcept {
    return Point2D(left.x / right, left.y / right);
}

std::istream& operator >> (std::istream& is, Point2D& p) {
   is >> p.x >> p.y;
   return is;
}

std::ostream& operator << (std::ostream& os, const Point2D& p) {
    os << p.x << ' ' << p.y;
    return os;
}

bool operator == (const Point2D& p1, const Point2D& p2) {
    if (boost::geometry::distance(p1, p2) < EPSILON) {
        return true;
    } else {
        return false;
    }
}

