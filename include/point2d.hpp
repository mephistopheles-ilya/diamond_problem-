#pragma once

#include <iostream>
#include <boost/geometry.hpp>

class Point2D {
public:
    double x = 0;
    double y = 0;
    int number_in_polyhedron = -1;
public:
    Point2D() noexcept = default;
    Point2D(double x, double y, int number_in_polyhedron = -1) noexcept: x(x), y(y), number_in_polyhedron(number_in_polyhedron) {}
    Point2D(const Point2D& other) noexcept: x(other.x), y(other.y), number_in_polyhedron(other.number_in_polyhedron) {}
    Point2D(Point2D&& other) noexcept: x(other.x), y(other.y), number_in_polyhedron(other.number_in_polyhedron) {}
    Point2D& operator = (const Point2D& other) & noexcept {
        x = other.x;
        y = other.y;
        number_in_polyhedron = other.number_in_polyhedron;
        return *this;
    }
    Point2D& operator = (Point2D&& other) & noexcept {
        x = other.x;
        y = other.y;
        number_in_polyhedron = other.number_in_polyhedron;
        return *this;
    }
    Point2D& operator += (const Point2D& other) & noexcept {
        x += other.x;
        y += other.y;
        return *this;
    }
    Point2D& operator -= (const Point2D& other) & noexcept {
        x -= other.x;
        y -= other.y;
        return *this;
    }
    Point2D& operator *= (double num) & noexcept {
        x *= num;
        y *= num;
        return *this;
    }
    Point2D& operator /= (double num) & noexcept {
        x /= num;
        y /= num;
        return *this;
    }
    friend bool operator == (const Point2D& p1, const Point2D& p2);
    ~Point2D() = default;

    friend Point2D operator + (const Point2D& left, const Point2D& right) noexcept;
    friend Point2D operator - (const Point2D& left, const Point2D& right) noexcept;
    friend Point2D operator * (const Point2D& left, double right) noexcept;
    friend Point2D operator * (double left, const Point2D& right) noexcept;
    friend Point2D operator / (const Point2D& left, double right) noexcept;
    friend std::istream& operator >> (std::istream& is, Point2D& p);
    friend std::ostream& operator << (std::ostream& os, const Point2D& p);
};


