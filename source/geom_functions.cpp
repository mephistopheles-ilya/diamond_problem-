#include "../include/geom_functions.hpp"
#include "../include/point2d.hpp"
#include "../include/point3d.hpp"
#include "../include/segment2d.hpp"
#include "../include/edge.hpp"

#include <boost/geometry/algorithms/correct.hpp>
#include <boost/geometry/algorithms/detail/distance/interface.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <ctime>
#include <iterator>
#include <random>

void uniform_grid_intersection(boost::geometry::model::linestring<Point2D>& line
        , boost::geometry::model::linestring<Point2D>& new_line
        , double dist_points) {
    double h = dist_points;
    Point2D point1 = line[0];
    Point2D point2 = line[1];
    for(size_t i = 1; i < line.size(); ++i) {
        point2 = line[i];
        Point2D vec = point2 - point1;
        double A = -vec.y;
        double B = vec.x;
        double C = point1.x * vec.y - point1.y * vec.x;
        if(std::abs(vec.x) < EPSILON) {
            new_line.push_back(point1);
            Point2D curent_point = point1;
            int y_index = (vec.y > 0) ? std::ceil(point1.y / h) : std::floor(point1.y / h);
            double next_y = y_index * h;
            Point2D next_point_y = Point2D(-(B * next_y + C) / A, next_y);

            bool y_flag = (vec.y > 0) ? ((next_y < point2.y) ? true : false) 
                : ((next_y > point2.y) ? true : false);
     
            while(y_flag) {
                new_line.push_back(next_point_y);
                curent_point = next_point_y;
                (vec.y > 0) ? (++y_index) : (--y_index);
                next_y = y_index * h;
                next_point_y = Point2D(-(B * next_y + C) / A, next_y);
                y_flag = (vec.y > 0) ? ((next_y < point2.y) ? true : false) 
                    : ((next_y > point2.y) ? true : false);
            }
        } else if (std::abs(vec.y) < EPSILON) {
            new_line.push_back(point1);
            Point2D curent_point = point1;
            int x_index = (vec.x > 0) ? std::ceil(point1.x / h) : std::floor(point1.x / h);
            double next_x = x_index * h;
            Point2D next_point_x = Point2D(next_x, -(A * next_x + C) / B);

            bool x_flag = (vec.x > 0) ? ((next_x < point2.x) ? true : false) 
                : ((next_x > point2.x) ? true : false);
     
            while(x_flag) {
                new_line.push_back(next_point_x);
                curent_point = next_point_x;
                (vec.x > 0) ? (++x_index) : (--x_index);
                next_x = x_index * h;
                next_point_x = Point2D(next_x, -(A * next_x + C) / B);
                x_flag = (vec.x > 0) ? ((next_x < point2.x) ? true : false) 
                    : ((next_x > point2.x) ? true : false);
            }
        } else {
            new_line.push_back(point1);
            Point2D curent_point = point1;
            int x_index = (vec.x > 0) ? std::ceil(point1.x / h) : std::floor(point1.x / h);
            int y_index = (vec.y > 0) ? std::ceil(point1.y / h) : std::floor(point1.y / h);
            double next_x = x_index * h;
            double next_y = y_index * h;
            Point2D next_point_x = Point2D(next_x, -(A * next_x + C) / B);
            Point2D next_point_y = Point2D(-(B * next_y + C) / A, next_y);

            bool x_flag = (vec.x > 0) ? ((next_x < point2.x) ? true : false) 
                : ((next_x > point2.x) ? true : false);
            bool y_flag = (vec.y > 0) ? ((next_y < point2.y) ? true : false) 
                : ((next_y > point2.y) ? true : false);
     
            while(x_flag || y_flag) {
                double dist_x = boost::geometry::distance(curent_point, next_point_x);
                double dist_y = boost::geometry::distance(curent_point, next_point_y);
                if(dist_x <=dist_y) {
                    new_line.push_back(next_point_x);
                    curent_point = next_point_x;
                    (vec.x > 0) ? (++x_index) : (--x_index);
                    next_x = x_index * h;
                    next_point_x = Point2D(next_x, -(A * next_x + C) / B);
                    x_flag = (vec.x > 0) ? ((next_x < point2.x) ? true : false) 
                        : ((next_x > point2.x) ? true : false);
                } else if (dist_y < dist_x) {
                    new_line.push_back(next_point_y);
                    curent_point = next_point_y;
                    (vec.y > 0) ? (++y_index) : (--y_index);
                    next_y = y_index * h;
                    next_point_y = Point2D(-(B * next_y + C) / A, next_y);
                    y_flag = (vec.y > 0) ? ((next_y < point2.y) ? true : false) 
                        : ((next_y > point2.y) ? true : false);
                }
#if 0
                else {
                    new_line.push_back(next_point_x);
                    curent_point = next_point_x;
                    ++x_index;
                    ++y_index;
                    (vec.x > 0) ? (++x_index) : (--x_index);
                    (vec.y > 0) ? (++y_index) : (--y_index);
                    next_point_x = Point2D(next_x, -(A * next_x + C) / B);
                    next_point_y = Point2D(-(B * next_y + C) / A, next_y);
                    x_flag = (vec.x > 0) ? ((next_x < point2.x) ? true : false) 
                        : ((next_x > point2.x) ? true : false);
                    y_flag = (vec.y > 0) ? ((next_y < point2.y) ? true : false) 
                        : ((next_y > point2.y) ? true : false);
                }
#endif
            }
        }
        point1 = point2;
    }
    new_line.push_back(point2);
}





Point2D rotate(Point2D p, double angle) {
    double a = cos(angle);
    double b = -sin(angle);
    double c = sin(angle);
    double d = cos(angle);
    return Point2D(p.x * a + p.y * b, p.x * c + p.y * d);
}

void create_small_shifts(boost::geometry::model::linestring<Point2D>& line, double sz_shift
        , std::string distribution, double sigma) {
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis_uni(-sz_shift, sz_shift);
    std::normal_distribution<> dis_normal(0, sigma);
    size_t end = line.size() - 1;
    if (distribution == "normal") {
        for(size_t i = 1; i < end; ++i) {
            Point2D direction = line[i] - line[i + 1];
            direction = Point2D(direction.y, -direction.x);
            double distance  = boost::geometry::distance(line[i], line[i + 1]);
            if (distance < EPSILON) {
                distance = boost::geometry::distance(line[i - 1], line[i + 2]);
                direction = line[i - 1] - line[i + 2];
            }
            direction = direction / distance;
            double l = dis_normal(gen);
            if (std::fabs(l) > sz_shift) {
                l = (l > 0) ? sz_shift : -sz_shift;
            }
            direction = direction * l;
            line[i] = line[i] + direction;
        }
    } else {
        for(size_t i = 1; i < end; ++i) {
            Point2D direction = line[i] - line[i + 1];
            direction = Point2D(direction.y, -direction.x);
            double distance  = boost::geometry::distance(line[i], line[i + 1]);
            if (distance < EPSILON) {
                distance = boost::geometry::distance(line[i - 1], line[i + 2]);
                direction = line[i - 1] - line[i + 2];
            }
            direction = direction / distance;
            direction = direction * dis_uni(gen);
            line[i] = line[i] + direction;
        }
    }

}


void spoil_and_get_protrusions(std::vector<Point2D>& v
        , std::vector<boost::geometry::model::polygon<Point2D, false, true, std::vector>>& vec_of_pol
        , std::vector<bool>& mask) {
    std::random_device rd;  
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis1(-boost::math::constants::pi<double>()/16
            , boost::math::constants::pi<double>()/16);
    std::uniform_real_distribution<> dis_small(__DISTANCE_BETWEEN_POINTS/2 
            , __DISTANCE_BETWEEN_POINTS);
    std::uniform_real_distribution<> dis_bigger(__DISTANCE_BETWEEN_POINTS
            , __DISTANCE_BETWEEN_POINTS * 2);

    boost::geometry::model::polygon<Point2D, false, true, std::vector> poly_of_spoiled_points;

    int flag = -1, step = 8;
    size_t j = 0, end = v.size() - 1;
    srand(time(NULL));

    for(size_t i = 1; i < end; ++i) {
        if ((i % 197) == 30 && (i + step < end) && (rand()%3 == 1)) {
            flag = (rand()%2) ? -1 : 1;
            if (flag == 1) mask.push_back(1);
            else mask.push_back(0);
            size_t index = i - 1;
            poly_of_spoiled_points.outer().push_back(v[i - 1]);
            for(j = i; j < (i + step); ++j) {
                Point2D direction = v[j] - v[j + 1];
                double distance = boost::geometry::distance(v[j], v[j + 1]);
                if (distance < EPSILON) {
                    distance = boost::geometry::distance(v[j - 1], v[j + 2]);
                    direction = v[j - 1] - v[j + 2];
                }
                direction = Point2D(direction.y, -direction.x) * flag;
                direction = direction / distance;
                //direction = rotate(direction, dis1(gen));
                Point2D vec = direction * (((j - i) < 3 || (j - i) > 5 ) 
                        ? dis_small(gen) : dis_bigger(gen));
                v[j] = v[j] + vec;
                poly_of_spoiled_points.outer().push_back(v[j]);
            }
            poly_of_spoiled_points.outer().push_back(v[j + 1]);
            poly_of_spoiled_points.outer().push_back(v[index]);
            vec_of_pol.push_back(poly_of_spoiled_points);
            boost::geometry::clear(poly_of_spoiled_points);
            i = j - 1;
       }
    }
}




void make_projection(double angle, std::vector<Point3D>& arr_points3d
        , boost::geometry::model::multi_point<Point2D, std::vector>& mpoint) {
    double sin_a = std::sin(angle);
    double cos_a = std::cos(angle);
    auto func = [sin_a, cos_a](Point3D& p) {
#if 0
        double length = std::sqrt(p.x * p.x + p.y * p.y);
        double sin_b = p.y / length;
        double cos_b = p.x / length;
        double cos_a_b = cos_a * cos_b + sin_a * sin_b;
        return Point2D(length * cos_a_b, p.z);
#endif
#if 0
        double length = std::sqrt(p.x * p.x + p.y * p.y);
        double sin_b = p.y / length;
        double cos_b = p.x / length;
        double sin_a_b = sin_a * cos_b - sin_b * cos_a;
        return Point2D(length * sin_a_b, p.z);
#endif
#if 1
        return Point2D(sin_a * p.x + cos_a * p.y, p.z, p.number_in_polyhedron);
#endif
    }; 
    std::ranges::transform(arr_points3d, std::back_inserter(mpoint), func);
}

Point2D find_minimal_y(boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull) {
    Point2D p_min_y = hull.outer()[0];
    for (auto& p : hull.outer()) {
        if (p.y < p_min_y.y) 
            p_min_y = p;
    }
    return p_min_y;
}

std::pair<double, double> find_borders_for_x(Point2D lowest
        , boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull) {
    double min_y_coord = lowest.y;
    double min_x_coord = 100;
    double max_x_coord = -100;
    for (auto& p : hull.outer()) {
        if ((p.y - min_y_coord) < BOTTOM_ERR) {
            if (p.x < min_x_coord)
                min_x_coord = p.x;
            if (p.x > max_x_coord)
                max_x_coord = p.x;
        }
    }
    return std::make_pair(min_x_coord, max_x_coord);
}

void projection_convex_hull(std::vector<Point3D>& arr_points3d
        , boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull, double angle) {
        boost::geometry::model::multi_point<Point2D, std::vector> mpoint;
        make_projection(angle, arr_points3d, mpoint);
        boost::geometry::convex_hull(mpoint, hull);
}

void projection_square_complexity(std::vector<Point3D>& arr_points3d
        , std::vector<Point3D>& arr_norm3d
        , std::vector<Edge>& arr_edges3d
        , boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull, double angle) {   
    std::vector<Segment2D> arr_segments2d;
    std::list<Segment2D> zero_proj2d;
    Point3D vec_proj(std::cos(angle), std::sin(angle), 0);
    for(auto& ed : arr_edges3d) {
        Point3D p1 = arr_points3d[ed.vert1_id];
        Point3D p2 = arr_points3d[ed.vert2_id];
        Point3D norm1 = arr_norm3d[ed.facet1_id];
        Point3D norm2 = arr_norm3d[ed.facet2_id];
        double dot_prod1 = vec_proj(norm1.norm());
        double dot_prod2 = vec_proj(norm2.norm());

        auto func = [sin_a = vec_proj.y, cos_a = vec_proj.x](Point3D& p) {
            double length = std::sqrt(p.x * p.x + p.y * p.y);
            double sin_b = p.y / length;
            double cos_b = p.x / length;
            double sin_a_b = sin_a * cos_b - sin_b * cos_a;
            return Point2D(length * sin_a_b, p.z);
        }; 
        if (dot_prod1 * dot_prod2 < 0) {
            Point2D proj1 = func(p1);
            Point2D proj2 = func(p2);
            arr_segments2d.push_back(Segment2D{proj1, proj2});
        }else if (dot_prod1 * dot_prod2 <= 0) {
            Point2D proj1 = func(p1);
            Point2D proj2 = func(p2);
            zero_proj2d.push_back(Segment2D{proj1, proj2});
        }
    }

    while(zero_proj2d.size()) {
        Segment2D seg = *zero_proj2d.begin();
        zero_proj2d.pop_front();
        for(auto it = zero_proj2d.begin(); it != zero_proj2d.end();) {
            if(seg.overlap(*it)) {
                seg = seg.seg_union(*it);
                zero_proj2d.erase(it);
                it = zero_proj2d.begin();
            } else {
                ++it;
            }
        }
        arr_segments2d.push_back(seg);
    }

    size_t all_points = arr_segments2d.size();
    Segment2D* current_seg = &arr_segments2d[0];
    hull.outer().push_back(current_seg -> p1);
    Point2D f_point = current_seg -> p1;
    Point2D start_point = current_seg -> p2;
    for(size_t i = 0; i < all_points; ++i) {
       for(auto& seg : arr_segments2d){
           if((boost::geometry::distance(seg.p1, start_point) < EPSILON) 
                   && (current_seg != &seg)) {
               start_point = seg.p2;
               hull.outer().push_back(seg.p1);
              if(boost::geometry::distance(seg.p1, f_point) < EPSILON && i > 0) {
                   i = all_points;
                   break;
               }
               current_seg = &seg;
               break;
           }else if((boost::geometry::distance(seg.p2, start_point) < EPSILON) 
                   && (current_seg != &seg)) {
               start_point = seg.p1;
               hull.outer().push_back(seg.p2);  
               if(boost::geometry::distance(seg.p2, f_point) < EPSILON && i > 0) {
                   i = all_points;
                   break;
               }
               current_seg = &seg;
               break;
           }
       }
    }

#if 1
    boost::geometry::model::polygon<Point2D, false, true, std::vector> tmp_hull;
    std::ranges::copy(hull.outer(), std::back_inserter(tmp_hull.outer()));
    boost::geometry::correct(tmp_hull);
    hull.clear();

    std::vector<Segment2D> segs;
    //TODO initial segment must not have any intersections
    Point2D first_point = tmp_hull.outer()[0];
    for(size_t i = 1; i < tmp_hull.outer().size(); ++i) {
        Segment2D seg = Segment2D(first_point, tmp_hull.outer()[i]);
        segs.push_back(seg);
        first_point = tmp_hull.outer()[i];
    }

    size_t i, j;
    for(i = 0; i < segs.size(); ++i) {
        for(j = i; j < segs.size(); ++j) {
            if (i != j && i != (j - 1) && i != (j + 1)) {
                auto [flag, point_of_intersect] = segs[i].intersect(segs[j]);
                if(flag && boost::geometry::distance(point_of_intersect, segs[i].p1) > EPSILON) {
                    hull.outer().push_back(segs[i].p1);
                    hull.outer().push_back(point_of_intersect);
                    hull.outer().push_back(segs[j].p2);
                    i = j;
                    break;
                }
            }
        }
        if(j == segs.size()) {
            hull.outer().push_back(segs[i].p1);
        }
    }
#endif



}



void hash_combine(size_t& seed, size_t hash_value) {
    seed ^= hash_value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}


void restore_polygon(const graph_t& graph, const Point2D& point, set_t& visited,
                     std::vector<Point2D>& polygon) {
    if (visited.find(point) == visited.end()) {
        visited.insert(point);
        polygon.push_back(point);
        for (const Point2D& q : graph.find(point)->second) {
            restore_polygon(graph, q, visited, polygon);
        }
    }
}

void restore_polygons(const std::vector<Segment2D>& segments,
                      std::vector<std::vector<Point2D>>& polygons) {
    graph_t graph;
    for (const Segment2D &s : segments) {
        graph[s.p1].push_back(s.p2);
        graph[s.p2].push_back(s.p1);
    }

    set_t visited;
    for (const auto& it : graph) {
        if (visited.find(it.first) == visited.end()) {
            polygons.push_back(std::vector<Point2D>());
            restore_polygon(graph, it.first, visited, polygons.back());
        }
    }
}


void projection_graph(std::vector<Point3D>& arr_points3d
        , std::vector<Point3D>& arr_norm3d
        , std::vector<Edge>& arr_edges3d
        , boost::geometry::model::polygon<Point2D, false, true, std::vector>& hull, double angle) {   
    std::vector<Segment2D> arr_segments2d;
    std::list<Segment2D> zero_proj2d;
    Point3D vec_proj(std::cos(angle), std::sin(angle), 0);
    for(auto& ed : arr_edges3d) {
        Point3D p1 = arr_points3d[ed.vert1_id];
        Point3D p2 = arr_points3d[ed.vert2_id];
        Point3D norm1 = arr_norm3d[ed.facet1_id];
        Point3D norm2 = arr_norm3d[ed.facet2_id];
        double dot_prod1 = vec_proj(norm1.norm());
        double dot_prod2 = vec_proj(norm2.norm());

        auto func = [sin_a = vec_proj.y, cos_a = vec_proj.x](Point3D& p) {
            double length = std::sqrt(p.x * p.x + p.y * p.y);
            double sin_b = p.y / length;
            double cos_b = p.x / length;
            double sin_a_b = sin_a * cos_b - sin_b * cos_a;
            return Point2D(length * sin_a_b, p.z);
        }; 

        if (dot_prod1 * dot_prod2 < 0) {
            Point2D proj1 = func(p1);
            Point2D proj2 = func(p2);
            arr_segments2d.push_back(Segment2D{proj1, proj2});
        }else if (dot_prod1 * dot_prod2 <= 0) {
            Point2D proj1 = func(p1);
            Point2D proj2 = func(p2);
            zero_proj2d.push_back(Segment2D{proj1, proj2});
        }
    }

    while(zero_proj2d.size()) {
        Segment2D seg = *zero_proj2d.begin();
        zero_proj2d.pop_front();
        for(auto it = zero_proj2d.begin(); it != zero_proj2d.end();) {
            if(seg.overlap(*it)) {
                seg = seg.seg_union(*it);
                zero_proj2d.erase(it);
                it = zero_proj2d.begin();
            } else {
                ++it;
            }
        }
        arr_segments2d.push_back(seg);
    }

    std::vector<std::vector<Point2D>> polygons;
    restore_polygons(arr_segments2d, polygons);
    size_t pol_index = 0;
    size_t max_sz = polygons[0].size();
    for(size_t i = 0; i < polygons.size(); ++i) {
        if(polygons[i].size() > max_sz) {
            max_sz = polygons[i].size();
            pol_index = i;
        }
    }
   // std::ranges::copy(polygons[0], std::back_inserter(hull.outer()));

#if 1
    boost::geometry::model::polygon<Point2D, false, true, std::vector> tmp_hull;
    std::ranges::copy(polygons[pol_index], std::back_inserter(tmp_hull.outer()));
    boost::geometry::correct(tmp_hull);

    std::vector<Segment2D> segs;
    //TODO initial segment must not have any intersections
    Point2D first_point = tmp_hull.outer()[0];
    for(size_t i = 1; i < tmp_hull.outer().size(); ++i) {
        Segment2D seg = Segment2D(first_point, tmp_hull.outer()[i]);
        segs.push_back(seg);
        first_point = tmp_hull.outer()[i];
    }

    size_t i, j;
    for(i = 0; i < segs.size(); ++i) {
        for(j = i; j < segs.size(); ++j) {
            if (i != j && i != (j - 1) && i != (j + 1)) {
                auto [flag, point_of_intersect] = segs[i].intersect(segs[j]);
                if(flag && boost::geometry::distance(point_of_intersect, segs[i].p1) > EPSILON) {
                    hull.outer().push_back(segs[i].p1);
                    hull.outer().push_back(point_of_intersect);
                    hull.outer().push_back(segs[j].p2);
                    i = j;
                    break;
                }
            }
        }
        if(j == segs.size()) {
            hull.outer().push_back(segs[i].p1);
        }
    }
#endif
}


