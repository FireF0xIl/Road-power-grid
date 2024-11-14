#ifndef GRAPH_H
#define GRAPH_H

class Map;

#include <set>
#include <utility>
#include "Config.h"
#include "Geometry.h"
#include "Map.h"
#include "GraphElement.h"


struct FakeVertex {
    std::set<size_t> next;
    std::set<size_t> real_next;
    Point point;
    size_t id;
    FakeVertex(Point point, size_t id) : point(point), id(id) {}
    FakeVertex(Point point) : point(point) {}

    bool operator<(const FakeVertex & p) const {
        return point < p.point;
    }
};

struct FakeGraph {
    std::map<FakeVertex, size_t> search;
    std::vector<FakeVertex> mp;

    size_t get_vertex(Point point) {
        FakeVertex v = FakeVertex(point);
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
    Graph() = default;
    std::vector<Vertex> vertex;
    std::set<Edge> edge;

    void build(const Map &map);

private:
    std::set<Point> point;
    std::map<Point, std::set<int>> point_polygon_map;
    std::map<std::pair<V_TYPE, Point>, size_t> point_vertex_id_map;

    void build_trench(const Map &map);
    void build_hdd(const Map &map);
    void check_hdd_connections(const FakeGraph &graph);
    void add_trench_point(const Map &map, int polygon_id, Point pt_from, Point pt_to);
    size_t get_vertex_id(V_TYPE v_type, const Point &e);
    void add_edge(V_TYPE v_type, EDGE_TYPE e_type, Edge &e);
};

#endif