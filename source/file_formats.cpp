#include <cstdint>
#include <stdio.h>
#include <tuple>
#include "../include/point2d.hpp"


bool save_con_file (const std::string& name, double pixel, std::tuple<double, double, double, double> plane 
        , const std::vector<std::vector<Point2D>>& contours, int ver) {
    std::string vers;
    switch (ver) {
        case 1: vers = "con1"; break;
        case 2: vers = "con2"; break;
        case 3: vers = "con3"; break;
        default: return false;
    }
    FILE* file = fopen(name.data(), "wb");
    if (fwrite(vers.data(), 4, 1, file) != 1) return false;
    if (fwrite(&pixel, 8, 1, file) != 1) return false;
    if (fwrite(&std::get<0>(plane), 8, 1, file) != 1 ) return false;
    if (fwrite(&std::get<1>(plane), 8, 1, file) != 1 ) return false;
    if (fwrite(&std::get<2>(plane), 8, 1, file) != 1 ) return false;
    if (fwrite(&std::get<3>(plane), 8, 1, file) != 1 ) return false;
    double scale = 0;
    unsigned int nc = static_cast<unsigned int>(contours.size());
    if(nc <= 1) {
        std::cerr << "Wrong amount of contours" << std::endl;
        return false;
    }
    if (ver == 3) {
        for (size_t i = 0; i < nc; ++i) {
            const std::vector<Point2D>& points = contours[i];
            unsigned int nv = static_cast<unsigned int>(points.size());
            for (size_t j = 0; j < nv; ++j) {
                scale = std::max(scale, std::fabs(points[j].x));
                scale = std::max(scale, std::fabs(points[j].y));
            }
        }
        if (scale <= 0) return false;
        scale /= 0x7ffffffd;
        if (fwrite(&scale, 8, 1, file) != 1 ) return false;
        scale = 1 / scale;
    }
    if (fwrite(&nc, 4, 1, file) != 1) return false;
    std::vector<float> buf;
    std::vector<int8_t> byte;
    for (size_t i = 0; i < nc; ++i) {
        const std::vector<Point2D>& points = contours[i];
        double angle = i * std::numbers::pi_v<double>/(nc - 1);
        if (fwrite(&angle, 8, 1, file) != 1 ) return false;
        unsigned int nv = points.size();
        if (fwrite(&nv, 4, 1, file) != 1 ) return false;
        switch (ver) {
            case 1:
                if (fwrite(points.data(), 16*nv, 1, file) != 1 ) return false;
                break;
            case 2:
                buf.resize(2 * nv);
                for (unsigned int j = 0; j < nv; ++j) {
                    buf[j + j] = static_cast<float>(points[j].x);
                    buf[j + j + 1] = static_cast<float>(points[j].y);
                }
                if (fwrite(buf.data(), 8*nv, 1, file) != 1 ) return false;
                break;
            case 3:
                byte.resize(6*nv);
                for (unsigned int j = 0; j < nv; ++j ) {
                    unsigned int k = 6 * j;
                    union IntByte {int32_t i; int8_t b[4]; } u;
                    u.i = static_cast<int>(scale * points[j].x);
                    byte[k] = u.b[1];
                    byte[k + 1] = u.b[2];
                    byte[k + 2] = u.b[3];
                    u.i = static_cast<int>(scale * points[j].y);
                    byte[k + 3] = u.b[1];
                    byte[k + 4] = u.b[2];
                    byte[k + 5] = u.b[3];
                }
                if (fwrite(byte.data(), 6*nv, 1, file) != 1 ) return false;
                break;
        }
    }
    return true;
}

// Функция возвращает номер версии файла или 0 в случае ошибки
int load_con_file (std::string name, double& pixel, std::tuple<double, double, double, double>& holder
        , std::vector<std::vector<Point2D>>& contours) {
    FILE* file = fopen(name.data(), "rb");
    char vers[4];
    if (fread(vers, 4, 1, file) != 1 ) return 0;
    if (vers[0] != 'c' || vers[1] != 'o' || vers[2] != 'n' ) return 0;
    switch (vers[3]) {
        case '1':
        case '2':
        case '3': break;
        default: return 0;
    }
    int ver = static_cast<int>(vers[3] - '0');
    if (fread(&pixel, 8, 1, file) != 1 ) return 0;
    if (fread(&std::get<0>(holder), 8, 1, file) != 1 ) return 0;
    if (fread(&std::get<1>(holder), 8, 1, file) != 1 ) return 0;
    if (fread(&std::get<2>(holder), 8, 1, file) != 1 ) return 0;
    if (fread(&std::get<3>(holder), 8, 1, file) != 1 ) return 0;
    double scale;
    if (ver == 3){
        if (fread(&scale, 8, 1, file) != 1 ) return false;
    }
    unsigned int nc;
    if (fread(&nc, 4, 1, file) != 1 ) return 0;
    contours.resize(nc);
    std::vector<float> buf;
    std::vector<int8_t> byte;
    for (unsigned int i = 0; i < nc; ++i) {
        std::vector<Point2D>& points = contours[i];
        double angle;
        if(fread(&angle, 8, 1, file) != 1 ) return 0;
        unsigned int nv;
        if (fread(&nv, 4, 1, file) != 1 ) return 0;
        points.resize(nv);
        switch (ver) {
            case 1:
                if (fread(points.data(), 16*nv, 1, file) != 1 ) return 0;
                break;
            case 2:
                float f[2];
                for (unsigned int j = 0; j < nv; ++j) {
                    if (fread(f, 8, 1, file) != 1 ) return 0;
                    points[j].x = static_cast<double>(f[0]);
                    points[j].y = static_cast<double>(f[1]);
                }
                break;
            case 3:
                    byte.resize(6*nv);
                    if(fread(byte.data(), 6*nv, 1, file) != 1 ) return false;
                    for (unsigned int j = 0; j < nv; ++j ) {
                        unsigned int k = 6 * j;
                        union IntByte { int32_t i; int8_t b[4]; } u;
                        u.b[0] = 0;
                        u.b[1] = byte[k];
                        u.b[2] = byte[k + 1];
                        u.b[3] = byte[k + 2];
                        points[j].x = scale * u.i;
                        u.b[1] = byte[k + 3];
                        u.b[2] = byte[k + 4];
                        u.b[3] = byte[k + 5];
                        points[j].y = scale * u.i;
                    }
        }
    }
    return ver;
}

