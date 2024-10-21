#pragma once

#include <iostream>

struct Edge {
    int vert1_id = 0, vert2_id = 0;
    int facet1_id = 0, facet2_id = 0;
    friend std::ostream& operator << (std::ostream& os, const Edge& ed);
};

