#include "points.h"

//function for Points on a plane
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


//functions for points in sapce
Point3D operator + (const Point3D& left, const Point3D& right) noexcept {
    return Point3D(left.x + right.x, left.y + right.y, left.z + right.z);
}

Point3D operator - (const Point3D& left, const Point3D& right) noexcept {
    return Point3D(left.x - right.x, left.y - right.y, left.z - right.z);
}

Point3D operator * (const Point3D& left, double right) noexcept {
    return Point3D(left.x * right, left.y * right, left.z * right);
}

Point3D operator * (double left, const Point3D& right) noexcept {
    return Point3D(right.x * left, right.y * left, right.z * left);
}

Point3D operator / (const Point3D& left, double right) noexcept {
    return Point3D(left.x / right, left.y / right, left.z / right);
}

std::istream& operator >> (std::istream& is, Point3D& p) {
   is >> p.x >> p.y >> p.z;
   return is;
}

std::ostream& operator << (std::ostream& os, const Point3D& p) {
    os << p.x << ' ' << p.y << ' ' << p.z;
    return os;
}


