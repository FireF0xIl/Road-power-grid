#include <iostream>
#include "Map.h"

void Map::fill_edges() {
    for (const auto &pol : polygons) {
        Point prev = pol.Points[0];
        for (size_t i = 1; i < pol.Points.size(); i++) {
            const Point cur = pol.Points[i];
            polygons_edges.emplace_back(prev, cur);
            prev = cur;
        }
    }
}

void Map::load_polygons(GeoJSONHandler &handler) {
    int polygon_id = 0;
    for (const auto &coordinates : handler.get_polygons_coordinates()) {
        Polygon cur(++polygon_id);
        for (auto &p: coordinates[0]) {
            const double x = p[0];
            const double y = p[1];
            cur.add_point(x, y);
        }
        polygons.emplace_back(cur);
    }
    fill_edges();
}

std::vector<std::pair<size_t, Point>> Map::find_segment_intersects(const Point p1, const Point p2) const {
    Point p;
    std::vector<std::pair<size_t, Point>> ans;
    for (size_t i = 0; i < polygons_edges.size(); i++) {
        const PolygonEdge edge = polygons_edges[i];
        if (segments_intersection(p1, p2, edge.point1, edge.point2, p)) {
            ans.emplace_back(i, p);
        }
    }
    return ans;
}


bool Map::check_edge_distance(const Point p1, const Point p2) const {
    std::vector<size_t> close_edges;
    const Point mid = mid_points(p1, p2);
    const double dist = p1.dist(p2);
    const double additional_dist = dist * 4.0;
    for (size_t i = 0; i < polygons_edges.size(); i++) {
        PolygonEdge cur_edge = polygons_edges[i];
        if (get_distance_from_segment(mid, cur_edge.point1, cur_edge.point2) < additional_dist) {
            close_edges.emplace_back(i);
        }
    }
    const int steps = static_cast<int>(ceil(dist / 10.0));
    for (int i = 1; i < steps; i++) {
        bool close = false;
        const Point cur_point = get_between_points(p1, p2, steps - i, i);
        if (!point_inside_polygons(cur_point, -1)) {
            for (size_t j : close_edges) {
                const PolygonEdge cur_edge = polygons_edges[j];
                if (get_distance_from_segment(cur_point, cur_edge.point1, cur_edge.point2) < MAX_DISTANCE_FROM_ROAD) {
                    close = true;
                    break;
                }
            }
            if (!close) {
                return false;
            }
        }
    }
    return true;
}

bool Map::point_inside_polygons(const Point p, int polygon_id) const {
    for (const auto &pp : polygons) {
        if (pp.id != polygon_id) {
            if (point_inside_polygon(p, pp)) {
                return true;
            }
        }
    }
    return false;
}


std::vector<Point> Map::get_intersection_points(int polygon_id, Point pt_from, Point pt_to) const {
    std::vector<Point>  pt_intersection;
    for (const auto& ps_outer: polygons) {
        if (polygon_id != ps_outer.id) {
            Point pt_outer_from = ps_outer.Points[ps_outer.Points.size() - 1];
            for (size_t i = 0; i < ps_outer.Points.size(); ++i) {
                Point pt_intersect;
                if (segments_intersection(pt_from, pt_to, pt_outer_from, ps_outer.Points[i], pt_intersect)) {
                    pt_intersection.push_back(pt_intersect);
                }
                pt_outer_from = ps_outer[i];
            }
        }
    }
    return pt_intersection;
}

double Map::get_alpha(Point vector, size_t id) const {
    return find_alpha(vector, getPerpendicularVector(polygons_edges[id].point1, polygons_edges[id].point2));
}