#include "../include/edge.h"

std::ostream& operator << (std::ostream& os, const Edge& ed) {
    os << ed.vert1_id << ' ' << ed.vert2_id << ' ' << ed.facet1_id << ' ' << ed.facet2_id;
    return os;
}

