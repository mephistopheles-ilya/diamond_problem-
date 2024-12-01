#!/bin/bash

#g++ -std=c++20 -O3 -DNDEBUG gorizontal_compare.cpp source/point2d.cpp -o gc
#g++ -std=c++20 -O3 -DNDEBUG gorizontal_compare_debug.cpp source/point2d.cpp source/segment2d.cpp -o gcd
#g++ -std=c++20 -O3 -DNDEBUG sym_diff.cpp source/point2d.cpp -o sd
g++ -std=c++20 -O3 -DNDEBUG -I ~/libs/boost_1_86_0 sym_diff.cpp source/point2d.cpp -o sd
#g++ -std=c++20 -O3 -DNDEBUG -I ~/libs/boost_1_86_0 sym_diff_debug.cpp source/point2d.cpp -o sddc
#g++ -std=c++20 -O3 -DNDEBUG splot.cpp source/point3d.cpp -lboost_program_options -lgmp
