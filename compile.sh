#!/bin/bash

g++ -std=c++20 -O3 -DNDEBUG gorizontal_compare.cpp source/point2d.cpp -o gc
g++ -std=c++20 -O3 -DNDEBUG sym_diff.cpp source/point2d.cpp -o sd
