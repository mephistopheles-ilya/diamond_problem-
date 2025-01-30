#include <iostream>
#include <tuple>
#include <vector>
#include <fstream>

#include "include/file_formats.hpp"
#include "include/point2d.hpp"

#include <boost/program_options.hpp>

int main(int argc, char* argv[]) {
    std::cerr << "Assuming that file contains all contours from 0 to 180 degrsees and they are soted from 0 to 180"
        << std::endl;
    std::string file;
    std::string directory;
    std::string type_of_file;
    boost::program_options::options_description desc("All options");
    desc.add_options()
        ("file", boost::program_options::value<std::string>(&file) 
         -> default_value("con_file"), "name of file in con format")
        ("directory", boost::program_options::value<std::string>(&directory)
         -> default_value("init_examples"), "name of dirictory whrere output files will be located")
        ("type_of_file", boost::program_options::value<std::string>(&type_of_file)
         -> default_value("angle"), "angle or number in postfix after Contour...") 
        ("help", "produce help message")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    double pixel = 0;
    std::tuple<double, double, double, double> holder;
    std::vector<std::vector<Point2D>> contours;
    int ver = 0;
    ver = load_con_file(file,pixel, holder, contours); 
    if (ver == 0) {
        std::cerr << "Cannot load con file" << std::endl;
        return 1;
    }
    size_t projections = contours.size() - 1;
    if (projections == 0) {
        std::cerr << "No contours in file" << std::endl;
        return 2;
    }
    for(size_t step = 0; step <= projections; ++step) {
        std::string file;
        if (type_of_file == "angle") {
            std::string angle_in_degrees = std::to_string(step * 180./projections);
            angle_in_degrees = angle_in_degrees.substr(0, angle_in_degrees.find(".") + 4);
            auto pos = angle_in_degrees.find(".");
            std::string zeros;
            while(pos < 3) {
                pos++;
                zeros += "0";
            }
            angle_in_degrees = zeros + angle_in_degrees;
            file = directory + std::string("/Contour") + angle_in_degrees + std::string(".txt");
        } else if (type_of_file == "number") {
            file = directory + std::string("/Contour") + std::to_string(step) + std::string(".txt");
        }
        std::ofstream print_res(file);
        if (!print_res.is_open()) {
            std::cerr << "problem with opening : " << file << std::endl;
        }
        std::ranges::copy(contours[step], std::ostream_iterator<Point2D>(print_res, "\n"));
    }
    return 0;
}



