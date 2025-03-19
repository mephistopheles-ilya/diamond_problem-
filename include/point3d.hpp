#pragma once

#include <iostream>
#include <cmath>

class Point3D {
public:
    double x = 0;
    double y = 0;
    double z = 0;
    int number_in_polyhedron = -1;
public:
    Point3D() noexcept = default;
    Point3D(double x, double y, double z, int number_in_polyhedron = -1) noexcept: x(x), y(y), z(z), number_in_polyhedron(number_in_polyhedron) {}
    Point3D(const Point3D& other) noexcept: x(other.x), y(other.y), z(other.z), number_in_polyhedron(other.number_in_polyhedron) {}
    Point3D(Point3D&& other) noexcept: x(other.x), y(other.y), z(other.z), number_in_polyhedron(other.number_in_polyhedron){}
    Point3D& operator = (const Point3D& other) & noexcept {
        x = other.x;
        y = other.y;
        z = other.z;
        number_in_polyhedron = other.number_in_polyhedron;
        return *this;
    }
    Point3D& operator = (Point3D&& other) & noexcept {
        x = other.x;
        y = other.y;
        z = other.z;
        number_in_polyhedron = other.number_in_polyhedron;
        return *this;
    }
    Point3D& operator += (const Point3D& other) & noexcept {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
    Point3D& operator -= (const Point3D& other) & noexcept {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }
    Point3D& operator *= (double num) & noexcept {
        x *= num;
        y *= num;
        z *= num;
        return *this;
    }
    Point3D& operator /= (double num) & noexcept {
        x /= num;
        y /= num;
        z /= num;
        return *this;
    }
    Point3D norm() const noexcept {
        double length = std::sqrt(x * x + y * y + z * z);
        return Point3D(x / length, y / length, z / length);
    }
    double operator() (const Point3D& other) const noexcept {
        return x * other.x + y * other.y + z * other.z;
    }
    ~Point3D() = default;

    friend Point3D operator + (const Point3D& left, const Point3D& right) noexcept;
    friend Point3D operator - (const Point3D& left, const Point3D& right) noexcept;
    friend Point3D operator * (const Point3D& left, double right) noexcept;
    friend Point3D operator * (double left, const Point3D& right) noexcept;
    friend Point3D operator / (const Point3D& left, double right) noexcept;
    friend bool operator == (const Point3D& left, const Point3D& right) noexcept;
    friend std::istream& operator >> (std::istream& is, Point3D& p);
    friend std::ostream& operator << (std::ostream& os, const Point3D& p);
};

