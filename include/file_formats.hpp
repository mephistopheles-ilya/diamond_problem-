#pragma once 

#include <tuple>
#include <vector>
#include "point2d.hpp"

bool save_con_file (const std::string& name, double pixel, std::tuple<double, double, double, double> plane 
        , const std::vector<std::vector<Point2D>>& contours, int ver); 
int load_con_file (std::string name, double& pixel, std::tuple<double, double, double, double>& holder
        , std::vector<std::vector<Point2D>>& contours); 

