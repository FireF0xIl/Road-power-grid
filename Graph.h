#ifndef GRAPH_H
#define GRAPH_H

#include <set>
#include <utility>
#include "nlohmann/json.hpp"
#include "Config.h"
#include "Geometry.h"

const double MIN_DIST = 10;

enum V_TYPE {
    TRENCH, HDD, TRENCH_TO_HDD, UNKNOWN
};

enum EDGE_TYPE {
    TRENCH_EDGE, HDD_EDGE, TRENCH_TO_HDD_EDGE, UNKNOWN_EDGE
};


struct GraphEdge {
    EDGE_TYPE o_type;
    size_t vertex_id;
    double cost;

    GraphEdge(EDGE_TYPE e_type, size_t vertex_id, double cost) : o_type(e_type), vertex_id(vertex_id), cost(cost) {}
};

class Vertex {
public:
    V_TYPE  o_type;
    Point position;
    std::vector<GraphEdge> graph_edges;

    Vertex(V_TYPE v_type, Point p) : o_type(v_type), position(p) {}

    bool operator<(const Vertex & p) const {
        if (o_type == p.o_type) {
            return position < p.position;
        } else {
            return o_type < p.o_type;
        }
    }

    V_TYPE type() const {
        return o_type;
    }

};

class Edge {
public:
    EDGE_TYPE o_type;
    Point point1, point2;
    double cost;

    Edge(EDGE_TYPE e_type, Point x, Point y): o_type(e_type), point1(x), point2(y) {
        cost = calculate_cost();
        if (point2 < point1) {
            std::swap(point1, point2);
        }
    }

    bool operator<(const Edge &edge) const {
        if (o_type != edge.o_type) {
            return o_type < edge.o_type;
        } else if (point1 == edge.point1) {
            return point2 < edge.point2;
        } else {
            return point1 < edge.point1;
        }
    }

    double getDist() const {
        return point1.dist(point2);
    }

    EDGE_TYPE type() const {
        return o_type;
    }

    double calculate_cost() const;

};

struct FakeVertex {
    std::set<size_t> next;
    std::set<size_t> real_next;
    Point point;
    V_TYPE type;
    size_t id;
    FakeVertex(Point &point, V_TYPE type, size_t id) : point(point), type(type), id(id) {}
    FakeVertex(Point &point, V_TYPE type) : point(point), type(type) {}

    bool operator<(const FakeVertex & p) const {
        if (type == p.type) {
            return point < p.point;
        } else {
            return type < p.type;
        }
    }
};

struct FakeGraph {
    std::map<FakeVertex, size_t> search;
    std::vector<FakeVertex> mp;

    size_t get_vertex(Point &point, V_TYPE type) {
        FakeVertex v = FakeVertex(point, type);
        size_t id;
        auto f = search.find(v);
        if (f == search.end()) {
            v.id = mp.size();
            id = v.id;
            search[v] = v.id;
            mp.emplace_back(v);
        } else {
            id = f->second;
        }
        return id;
    }

    void add_edge(size_t id1, size_t id2) {
        mp[id1].next.insert(id2);
        mp[id2].next.insert(id1);
    }
};

class Graph {
public:
    Graph();
    std::vector<Polygon> polygons;
    std::vector<Vertex> vertex;
    std::set<Edge> edge;

    void load();
    void build();

    void save();
private:
    nlohmann::json jf;

    std::set<Point> point;
    std::map<Point, std::set<int>> point_polygon_map;
    std::map<std::pair<V_TYPE, Point>, size_t> point_vertex_id_map;
    std::map<Vertex, size_t> vertex_id_map;

    void build_trench();
    void build_hdd();
    void check_hdd_connections(const FakeGraph &graph);
    bool point_inside_polygons(const Point &p, int polygon_id = -1);
    bool hdd_roads_intersection(const Point &p1, const Point &p2);
    bool segment_polygons_intersection(const Point &p1, const Point &p2);
    bool point_inside_polygons(const Point &p, const std::set<int> &polygon_ids = std::set<int>());
    std::vector<Point> get_intersection_points(int polygon_id, Point pt_from, Point pt_to);
    void add_trench_point(int polygon_id, Point &pt_from, const Point &pt_to);
    size_t get_vertex_id(V_TYPE v_type, const Point &e);
    void add_edge(V_TYPE v_type, EDGE_TYPE e_type, Edge &e);
    void build_additional_trenches(size_t vertex_start);
    bool check_edge_distance(const std::vector<Edge> &polygons_edges, const Point &p1, const Point &p2);
};

#endif