#pragma once

#include <iostream>

class Point2D {
public:
    double x = 0;
    double y = 0;
public:
    Point2D() noexcept = default;
    Point2D(double x, double y) noexcept: x(x), y(y) {}
    Point2D(const Point2D& other) noexcept: x(other.x), y(other.y) {}
    Point2D(Point2D&& other) noexcept: x(other.x), y(other.y) {}
    Point2D& operator = (const Point2D& other) & noexcept {
        x = other.x;
        y = other.y;
        return *this;
    }
    Point2D& operator = (Point2D&& other) & noexcept {
        x = other.x;
        y = other.y;
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
    ~Point2D() = default;

    friend Point2D operator + (const Point2D& left, const Point2D& right) noexcept;
    friend Point2D operator - (const Point2D& left, const Point2D& right) noexcept;
    friend Point2D operator * (const Point2D& left, double right) noexcept;
    friend Point2D operator * (double left, const Point2D& right) noexcept;
    friend Point2D operator / (const Point2D& left, double right) noexcept;
    friend std::istream& operator >> (std::istream& is, Point2D& p);
    friend std::ostream& operator << (std::ostream& os, const Point2D& p);
};



class Point3D {
public:
    double x = 0;
    double y = 0;
    double z = 0;
public:
    Point3D() noexcept = default;
    Point3D(double x, double y, double z) noexcept: x(x), y(y), z(z) {}
    Point3D(const Point3D& other) noexcept: x(other.x), y(other.y), z(other.z) {}
    Point3D(Point3D&& other) noexcept: x(other.x), y(other.y), z(other.z) {}
    Point3D& operator = (const Point3D& other) & noexcept {
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }
    Point3D& operator = (Point3D&& other) & noexcept {
        x = other.x;
        y = other.y;
        z = other.z;
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
    ~Point3D() = default;

    friend Point3D operator + (const Point3D& left, const Point3D& right) noexcept;
    friend Point3D operator - (const Point3D& left, const Point3D& right) noexcept;
    friend Point3D operator * (const Point3D& left, double right) noexcept;
    friend Point3D operator * (double left, const Point3D& right) noexcept;
    friend Point3D operator / (const Point3D& left, double right) noexcept;
    friend std::istream& operator >> (std::istream& is, Point3D& p);
    friend std::ostream& operator << (std::ostream& os, const Point3D& p);

};

