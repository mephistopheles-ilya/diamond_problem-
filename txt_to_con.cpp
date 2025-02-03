#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <tuple>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <cstdio>
#include "include/geom_functions.hpp"

#include "include/file_formats.hpp"
#include "include/point2d.hpp"

#include <boost/program_options.hpp>

int main(int argc, char* argv[]) {
    std::cerr << "Assuming that directory contains all contours from 0 to 180 degrees" << std::endl;
    std::string name_of_con_file;
    std::string directory;
    std::string type_of_file;
    std::string model_in_txt;
    boost::program_options::options_description desc("All options");
    desc.add_options()
        ("name_of_con_file", boost::program_options::value<std::string>(&name_of_con_file) 
         -> default_value("con_file"), "name of file in con format")
        ("directory", boost::program_options::value<std::string>(&directory)
         -> default_value("./init_examples"), "name of dirictory whrere files with points are located (full path)")
        ("type_of_file", boost::program_options::value<std::string>(&type_of_file)
         -> default_value("con1"), "con1 or con2 or con3") 
        ("model_in_txt", boost::program_options::value<std::string>(&model_in_txt) 
         -> default_value("InitialModels/InitialModel5"), "3d model of contours to get holder")
        ("help", "produce help message")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    std::vector<std::vector<Point2D>> contours;

    std::vector<std::string> filenames;
    std::filesystem::path p(directory);
    if (std::filesystem::exists(p) && std::filesystem::is_directory(p)) {
        for (const auto& entry : std::filesystem::directory_iterator(p)) {
            filenames.push_back(entry.path().filename());
        }
    } else {
        std::cerr << "Directory not found\n";
        return 1;
    }

    std::sort(filenames.begin(), filenames.end());
    unsigned int nc = 0;
    for(auto& file: filenames) {
        if(file.find("Contour") != std::string::npos) {
            ++nc;
        }
    }
    if (nc == 0) {
        std::cerr << "No files in deirectory" << std::endl;
        return 2;
    }
    size_t step = 0;
    contours.resize(nc);
    for(auto& file: filenames) {
        if(file.find("Contour") != std::string::npos) {
            std::ifstream istr(directory + std::string("/") + file);
            std::istream_iterator<Point2D> it(istr);
            std::copy(it, std::istream_iterator<Point2D>{}, std::back_inserter(contours[step]));
            ++step;
        }
    }

    double pixel = 0;
    size_t amount_of_points = 0;
    for(size_t i = 0; i < nc; ++i) {
        std::vector<Point2D>& points = contours[i];
        unsigned int nv = points.size();
        amount_of_points += nv;
        for(size_t j = 0; j < (nv - 1); ++j) {
            pixel += boost::geometry::distance(points[i], points[i + 1]);
        }
    }
    pixel /= amount_of_points;
    pixel = std::sqrt(pixel);
    int ver = 0;
    if (type_of_file == "con1") {
        ver = 1;
    }
    if (type_of_file == "con2") {
        ver = 2;
    }
    if (type_of_file == "con3") {
        ver = 3;
    }
    std::ifstream iss(model_in_txt);
    double min_z = std::numeric_limits<double>().max();
    if (iss.is_open() == false) {
        min_z = -1;
    } else {
        std::string line;
        std::string start("vertices:");
        std::string stop("#");
        std::getline(iss, line);
        while(line.find(start) == std::string::npos) {
            std::getline(iss, line);
        }
        std::getline(iss, line);
        std::getline(iss, line);
        while(line.find(stop) == std::string::npos) {
            int number = 0;
            double x = 0, y = 0, z = 0;
            sscanf(line.data(), "%d %lf %lf %lf", &number, &x, &y, &z);
            min_z = std::min(z, min_z);
            std::getline(iss, line);
        }
    }

    bool ret;
    ret = save_con_file(name_of_con_file, pixel, std::make_tuple(0, 0, 1, -min_z), contours, ver);
    if (ret == false) {
        std::cerr << "Cannot write con file" << std::endl;
    }


}

