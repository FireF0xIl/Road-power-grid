#ifndef HDD_GRAPHELEMENT_H
#define HDD_GRAPHELEMENT_H

#include "Geometry.h"
#include "Config.h"

enum class V_TYPE {
    TRENCH, HDD, UNKNOWN
};

enum class EDGE_TYPE {
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

    auto operator<=>(const Vertex &p) const {
        if (o_type == p.o_type) {
            return position <=> p.position;
        } else if (o_type < p.o_type) {
            return std::strong_ordering::less;
        } else {
            return std::strong_ordering::greater;
        }
    }

    bool operator ==(const Vertex &p) const {
        return (*this <=> p) == std::strong_ordering::equivalent;
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

    auto operator<=>(const Edge &edge) const {
        if (o_type == edge.o_type) {
            if (auto cmp = point1 <=> edge.point1; cmp != std::strong_ordering::equivalent)
                return cmp;
            return point2 <=> edge.point2;
        } else if (o_type < edge.o_type) {
            return std::strong_ordering::less;
        } else {
            return std::strong_ordering::greater;
        }
    }

    bool operator ==(const Edge &edge) const {
        return (*this <=> edge) == std::strong_ordering::equivalent;
    }

    double get_dist() const {
        return point1.dist(point2);
    }

    EDGE_TYPE type() const {
        return o_type;
    }

    double calculate_cost() const {
        double value;
        if (o_type == EDGE_TYPE::TRENCH_EDGE) {
            value = get_dist() * Config::trench_cost;
        } else if(o_type == EDGE_TYPE::HDD_EDGE) {
            value = get_dist() * Config::hdd_cost;
        } else if (o_type == EDGE_TYPE::TRENCH_TO_HDD_EDGE) {
            value = Config::trench_hdd_cost;
        } else {
            value = 0.0;
        }
        return value;
    }
};

#endif
