#pragma once

#include "point2d.h"

struct Segment2D {
    Point2D p1, p2;
    bool check_line(const Segment2D& other);
    bool overlap(const Segment2D& other);
    Segment2D seg_union(const Segment2D& other);
    friend std::ostream& operator << (std::ostream& os, const Segment2D& seg);
};

