#ifndef HDD_MAP_H
#define HDD_MAP_H

class GeoJSONHandler;

#include <vector>
#include "GeoJSONHandler.h"
#include "GraphElement.h"

struct PolygonEdge {
public:
    Point point1, point2;

    PolygonEdge(Point x, Point y): point1(x), point2(y) {
        if (point2 < point1) {
            std::swap(point1, point2);
        }
    }
};

class Map {
public:
    Map() = default;
    std::vector<Polygon> polygons;
    std::vector<PolygonEdge> polygons_edges;

    double get_alpha(Point vector, size_t id) const;
    void load_polygons(GeoJSONHandler &GeoJSONHandler);
    std::vector<std::pair<size_t, Point>> find_segment_intersects(Point p1, Point p2) const;

    bool point_inside_polygons(Point p, int polygon_id) const;
    bool check_edge_distance(Point p1, Point p2) const;
    std::vector<Point> get_intersection_points(int polygon_id, Point pt_from, Point pt_to) const;

private:
    void fill_edges();
};


#endif
