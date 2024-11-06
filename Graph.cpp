#include <iostream>
//#include "Geometry.h"
#include "Graph.h"
#include "GeoJSONUtility.h"
#include <stack>
#include <utility>

static Config config;

Graph::Graph() {
    config.read();
}

void Graph::load() {
    if (config.input.empty()) {
        return;
    }
    jf = GeoJSONReader(config.input);

    int polygon_id = 0;
    for (auto &r: jf["features"]) {
        for (auto &coordinates: r["geometry"]["coordinates"]) {
            Polygon cur(++polygon_id);
            for (auto &p: coordinates) {
                double x = p[0];
                double y = p[1];
                cur.add_point(x, y);
            }
            polygons.emplace_back(cur);
        }
    }
}

void Graph::build() {
    build_trench();
    build_hdd();
}

void Graph::build_trench() {
    std::cout << "Trench build start" << std::endl;
    for (const auto& ps: polygons) {
        Point pt_from = ps.Points[ps.Points.size() - 1];

        if (!point_inside_polygons(pt_from, ps.id)) {
            auto result = point.insert(pt_from);
            if (result.second) {
                point_polygon_map.insert({pt_from, {ps.id}});
            } else {
                auto value = point_polygon_map[pt_from];
                value.insert(ps.id);
                point_polygon_map[pt_from] = value;
            }
        }

        for (size_t i = 0; i < ps.Points.size(); ++i) {
            Point pt_to = ps.Points[i];
            std::vector<Point> pt_intersect = get_intersection_points(ps.id, pt_from, pt_to);
            std::sort(pt_intersect.begin(), pt_intersect.end());
            if (!pt_intersect.empty()) {
                if (ps.Points[i] < pt_from) {
                    std::reverse(pt_intersect.begin(), pt_intersect.end());
                }
                for (const auto pt_i: pt_intersect) {
                    add_trench_point(ps.id, pt_from, pt_i);
                    pt_from = pt_i;
                }
                if (!point_inside_polygons(pt_to, ps.id)) {
                    add_trench_point(ps.id, pt_from, pt_to);
                } else {
                }
            } else {
                if (!point_inside_polygons(pt_to, ps.id)) {
                    add_trench_point(ps.id, pt_from, pt_to);
                }
            }
            pt_from = ps.Points[i];
        }
    }
    for (auto e: edge) {
        add_edge(TRENCH, TRENCH_EDGE, e);
    }

    std::cout << "Trench build finish" << std::endl;
}

void Graph::add_edge(V_TYPE v_type, EDGE_TYPE e_type, Edge &e) {
    size_t v1_id = get_vertex_id(v_type, e.point1);
    size_t v2_id = get_vertex_id(v_type, e.point2);
    vertex[v1_id].graph_edges.emplace_back(e_type, v2_id, e.cost);
    vertex[v2_id].graph_edges.emplace_back(e_type, v1_id, e.cost);
}

size_t Graph::get_vertex_id(V_TYPE v_type, const Point &p) {
    size_t v_id;
    if (point_vertex_id_map.find({v_type, p}) != point_vertex_id_map.end()) {
        v_id = point_vertex_id_map[{v_type, p}];
    } else {
        v_id = vertex.size();
        Vertex vert = Vertex(v_type, p);
        vertex.emplace_back(vert);
        point_vertex_id_map[{v_type, p}] = v_id;
        vertex_id_map[vert] = v_id;
    }
    return v_id;
}

void Graph::add_trench_point(int polygon_id, Point &pt_from, const Point &pt_to) {
    if (point_inside_polygons(Point({(pt_from.x + pt_to.x) / 2.0, (pt_from.y + pt_to.y) / 2.0}), polygon_id)) {
        auto result = point.insert(pt_to);
        if (result.second) {
            point_polygon_map.insert({pt_to, {polygon_id}});
        } else {
            auto value = point_polygon_map[pt_to];
            value.insert(polygon_id);
            point_polygon_map[pt_to] = value;
        }
    } else {
        auto pt = splice_segment(pt_from, pt_to);
        Point pt_cur = pt_from;
        for (size_t i = 1; i < pt.size(); ++i) {
            auto result = point.insert(pt[i]);
            if (result.second) {
                point_polygon_map.insert({pt[i], {polygon_id}});
            } else {
                auto value = point_polygon_map[pt[i]];
                value.insert(polygon_id);
                point_polygon_map[pt[i]] = value;
            }

            edge.insert(Edge(TRENCH_EDGE, pt_cur, pt[i]));
            pt_cur = pt[i];
        }
    }
}

// (-c*dy, c*dx) - perpendicular, c - coefficient
std::pair<Rectangle, Rectangle> build_rectangle(const Point &point1, const Point &point2, size_t id1, size_t id2) {
    double dx = point1.x - point2.x;
    double dy = point1.y - point2.y;
    // norm

    double base = hypot(dx, dy);
    dx = dx / base * MIN_DIST;
    dy = dy / base * MIN_DIST;
    Rectangle first = {
            point1,
            point2,
            Point(point1.x - dy, point1.y + dx),
            Point(point2.x - dy, point2.y + dx),
            id1,
            id2
    };
    Rectangle second = {
            point1,
            point2,
            Point(point1.x + dy, point1.y - dx),
            Point(point2.x + dy, point2.y - dx),
            id1,
            id2
    };

    return std::make_pair(first, second);
}

void Graph::check_hdd_connections(const FakeGraph &graph) {
    std::vector<int> visited(graph.mp.size(), -1);
    std::map<int, int> col;
    int colour = 0;
    int max = std::numeric_limits<int>::min();
    int chosen_colour;
    for (size_t i = 0; i < graph.mp.size(); i++) {
        if (graph.mp[i].type == HDD && visited[i] == -1) {
            int cur_colour = colour++;
            std::stack<size_t> stack;
            stack.push(i);
            int count = 0;
            while (!stack.empty()) {
                size_t cur_i = stack.top();
                stack.pop();
                visited[cur_i] = cur_colour;
                count++;
                for (auto j : graph.mp[cur_i].next) {
                    if (graph.mp[j].type == HDD && visited[j] == -1) {
                        stack.push(j);
                    }
                }
            }
            if (max < count) {
                max = count;
                chosen_colour = cur_colour;
            }
            col[cur_colour] = count;
        }
    }
    std::vector<size_t> hdd;
    for (size_t i = 0; i < visited.size(); i++) {
        if (visited[i] == chosen_colour) {
            hdd.emplace_back(i);
        }
    }
    std::vector<size_t> trench;
    for (auto i : hdd) {
        for (auto e : graph.mp[i].next) {
            if (graph.mp[e].type == TRENCH) {
                trench.emplace_back(e);
            }
        }
    }

    std::map<size_t, size_t> translation; // fake to real
    for (auto i : hdd) {
        translation[graph.mp[i].id] = get_vertex_id(HDD, graph.mp[i].point);
    }
    for (auto i : trench) {
        translation[graph.mp[i].id] = get_vertex_id(TRENCH, graph.mp[i].point);
    }
    for (size_t i = 0; i < hdd.size(); i++) {
        size_t id = translation[hdd[i]];
        for (auto e : graph.mp[hdd[i]].next) {
            double dist = graph.mp[e].point.dist(graph.mp[hdd[i]].point);
            double cost;
            EDGE_TYPE type;
            if (graph.mp[e].type == HDD) {
                type = HDD_EDGE;
                cost = dist * config.hdd_cost;
            } else {
                type = TRENCH_TO_HDD_EDGE;
                cost = config.trench_hdd_cost;
            }
            vertex[id].graph_edges.emplace_back(type, translation[e], cost);
        }
        for (auto e : graph.mp[hdd[i]].real_next) {
            EDGE_TYPE type = TRENCH_TO_HDD_EDGE;
            double cost = config.trench_hdd_cost;
            vertex[e].graph_edges.emplace_back(type, id, cost);
            vertex[id].graph_edges.emplace_back(type, e, cost);
        }
    }
    for (size_t i = 0; i < trench.size(); i++) {
        size_t id = translation[trench[i]];
        for (auto e : graph.mp[trench[i]].next) {
            vertex[id].graph_edges.emplace_back(TRENCH_TO_HDD_EDGE, translation[e], config.trench_hdd_cost);
        }
        for (auto e : graph.mp[trench[i]].real_next) {
            EDGE_TYPE type = TRENCH_EDGE;
            double cost = vertex[e].position.dist(vertex[id].position) * config.trench_cost;
            vertex[e].graph_edges.emplace_back(type, id, cost);
            vertex[id].graph_edges.emplace_back(type, e, cost);
        }
    }

    for (size_t i = 0; i < vertex.size(); i++) {
        auto q = vertex[i];
        if (vertex[i].o_type == HDD) {
            for (auto &e : vertex[i].graph_edges) {
                if (e.vertex_id > i && e.o_type == HDD_EDGE) {
                    edge.insert(Edge(HDD_EDGE, vertex[i].position, vertex[e.vertex_id].position));
                }
            }
        }
    }

    for (size_t i = 0; i < vertex.size(); i++) {
        auto q = vertex[i];
        if (vertex[i].o_type == TRENCH) {
            for (auto &e : vertex[i].graph_edges) {
                if (e.vertex_id > i && e.o_type == TRENCH_EDGE) {
                    edge.insert(Edge(TRENCH_EDGE, vertex[i].position, vertex[e.vertex_id].position));
                }
            }
        }
    }
    for (size_t i = 0; i < vertex.size(); i++) {
        auto q = vertex[i];
        if (vertex[i].o_type == TRENCH || vertex[i].o_type == HDD) {
            for (auto &e : vertex[i].graph_edges) {
                if (e.vertex_id > i && e.o_type == TRENCH_TO_HDD_EDGE) {
                    edge.insert(Edge(TRENCH_TO_HDD_EDGE, vertex[i].position, vertex[e.vertex_id].position));
                }
            }
        }
    }
}

void Graph::build_additional_trenches(const size_t vertex_start) {
    for (size_t i = 0; i < vertex.size(); i++) {
        if (vertex[i].o_type == TRENCH) {
            auto v1 = vertex[i];
            for (size_t j = vertex_start; j < vertex.size(); j++) {
                if (vertex[j].o_type == TRENCH) {
                    auto v2 = vertex[j];
                    double dist = v1.position.dist(v2.position);
                    if (dist < 1.5 * DENSITY) {
                        Edge trench_edge = Edge(TRENCH_EDGE, v1.position, v2.position);
                        if (edge.find(trench_edge) != edge.end()) {
                            continue;
                        }
                        if (!segment_polygons_intersection(v1.position, v2.position)) {
                            double cost = dist * config.trench_cost;
                            v1.graph_edges.emplace_back(TRENCH_EDGE, j, cost);
                            v2.graph_edges.emplace_back(TRENCH_EDGE, i, cost);
                            edge.insert(trench_edge);
                        }
                    }
                }
            }
        }
    }
}

void write_vertex(FakeGraph &local_graph,
                  Point &point,
                  int k1,
                  int k2,
                  size_t h_id,
                  size_t road1,
                  size_t road2) {
    if (k1 == 0) {
        // на дороге
        if (k2 == 0) {
            local_graph.mp[h_id].real_next.insert(road1);
        } else {
            local_graph.mp[h_id].real_next.insert(road2);
        }
    } else {
        // не на дороге
        auto t_id = local_graph.get_vertex(point, TRENCH);
        local_graph.add_edge(t_id, h_id);
        if (k2 == 0) {
            local_graph.mp[t_id].real_next.insert(road1);
        } else {
            local_graph.mp[t_id].real_next.insert(road2);
        }
    }
}

void add_hdd_edge(FakeGraph &local_graph,
                  Point &p1,
                  Point &p2,
                  const Rectangle &r1,
                  const Rectangle &r2,
                  int k1,
                  int k2,
                  int k3,
                  int k4) {
    auto h1_id = local_graph.get_vertex(p1, HDD);
    auto h2_id = local_graph.get_vertex(p2, HDD);
    local_graph.add_edge(h1_id, h2_id);

    size_t road1 = r1.id1;
    size_t road2 = r1.id2;
    size_t road3 = r2.id1;
    size_t road4 = r2.id2;
    write_vertex(local_graph, p1, k1, k2, h1_id, road1, road2);
    write_vertex(local_graph, p2, k3, k4, h2_id, road3, road4);
}

void Graph::build_hdd() {
    std::cout << "HDD build start" << std::endl;
    std::vector<Rectangle> rectangles_new;
    // считаем, что отступ от дороги достаточно мал, чтобы не попасть на другую дорогу
    std::vector<int> visited(vertex.size(), -1);

    for (size_t i = 0; i < vertex.size(); i++) {
        if (visited[i] == -1) {
            visited[i] = 0;
            Vertex cur = vertex[i];
            for (const auto &j : cur.graph_edges) {
                if (visited[j.vertex_id] == -1) {
                    continue;
                }
                if (cur.position == vertex[j.vertex_id].position) {
                    continue;
                }
                auto pair = build_rectangle(cur.position, vertex[j.vertex_id].position, i, j.vertex_id);
                if (!point_inside_polygons(pair.first.mid_point(), -1)) {
                    rectangles_new.emplace_back(pair.first);
                }
                if (!point_inside_polygons(pair.second.mid_point(), -1)) {
                    rectangles_new.emplace_back(pair.second);
                }
            }
        }
    }

    std::set<Vertex> v_local;
    FakeGraph local_graph = FakeGraph();
    for (int i = 0; i < rectangles_new.size(); i++) {
        auto r1 = rectangles_new[i];
        // не перебирать прямоугольники, а искать пересечение перпендикуляров к дороге
        for (int j = i + 1; j < rectangles_new.size(); j++) {
            auto r2 = rectangles_new[j];
            double rect_dist = r1.center().dist(r2.center());
            if (rect_dist <= config.min_hdd_distance - DENSITY || config.max_hdd_distance + DENSITY <= rect_dist) {
                continue;
            }
            for (int k1 = 0; k1 <= 1; k1++) {
                Point pp1 = get_between_points(r1.base1, r1.perp1, 1 - k1, k1);
                Point pp2 = get_between_points(r1.base2, r1.perp2, 1 - k1, k1);
                // рассматривается только самая дальняя точка
                for (int k2 = 0; k2 <= 1; k2++) {
                    Point p1 = get_between_points(pp1, pp2, 1 - k2, k2);
                    for (int k3 = 0; k3 <= 1; k3++) {
                        Point pp3 = get_between_points(r2.base1, r2.perp1, 1 - k3, k3);
                        Point pp4 = get_between_points(r2.base2, r2.perp2, 1 - k3, k3);
                        for (int k4 = 0; k4 <= 1; k4++) {
                            Point p2 = get_between_points(pp3, pp4, 1 - k4, k4);
                            double dist = p1.dist(p2);
                            if (config.min_hdd_distance <= dist && dist <= config.max_hdd_distance) {
                                Point t;
                                if (segments_intersection(p1, p2, r1.base1, r1.base2, t)
                                    && segments_intersection(p1, p2, r2.base1, r2.base2, t)) {
                                    Point p1p2 = get_vector(p1, p2);
                                    double alpha1 = find_alpha(p1p2, r1.get_perpendicular_vector());
                                    double alpha2 = find_alpha(p1p2, r2.get_perpendicular_vector());
                                    // Построение по дороге ИЛИ посроение вне дорог
                                    if (alpha1 < config.alpha
                                        && alpha2 < config.alpha
                                        && !point_inside_polygons(p1, point_polygon_map[p1])
                                        && !point_inside_polygons(p2, point_polygon_map[p2])) {
                                        add_hdd_edge(local_graph, p1, p2, r1, r2, k1, k2, k3, k4);
                                    }
                                } else if (!hdd_roads_intersection(p1, p2)) {
                                    add_hdd_edge(local_graph, p1, p2, r1, r2, k1, k2, k3, k4);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    size_t vertex_size = vertex.size();
    check_hdd_connections(local_graph);
    build_additional_trenches(vertex_size);
    std::cout << "HDD build finish" << std::endl;
}

void Graph::save() {
    GeoJSONWriter(jf, vertex, TRENCH, "trench_");
    GeoJSONWriter(jf, vertex, HDD, "hdd_");
    GeoJSONWriter(jf, edge, TRENCH_EDGE, "trench_");
    GeoJSONWriter(jf, edge, HDD_EDGE, "hdd_");
    GeoJSONWriter(jf, edge, TRENCH_TO_HDD_EDGE, "trench_to_hdd_");

}


bool segments_intersection(const Point &p1, const Point &q1,
                           const Edge &edge,
                           Point &intersection) {
    return segments_intersection(p1, q1, edge.point1, edge.point2, intersection);
}

bool Graph::point_inside_polygons(const Point &p, int polygon_id) {
    for (const auto &pp : polygons) {
        if (pp.id != polygon_id) {
            if (point_inside_polygon(p, pp)) {
                return true;
            }
        }
    }
    return false;
}

bool Graph::point_inside_polygons(const Point &p, const std::set<int> &polygon_ids) {
    for (const auto &pp : polygons) {
        if (polygon_ids.count(pp.id) == 0) {
            if (point_inside_polygon(p, pp)) {
                return true;
            }
        }
    }
    return false;
}

std::pair<Point, Point> get_close_points(const Point &p1_, const Point &p2_) {
    double dx = p2_.x - p1_.x;
    double dy = p2_.y - p1_.y;
    // norm
    double base = hypot(dx, dy);
    dx = dx / base;
    dy = dy / base;
    Point p1 = {p1_.x + dx, p1_.y + dy};
    Point p2 = {p2_.x - dx,  p2_.y - dy};
    return std::make_pair(p1, p2);
}

bool Graph::hdd_roads_intersection(const Point &p1_, const Point &p2_) {
    Point p;
    auto pair = get_close_points(p1_, p2_);
    Point p1 = pair.first;
    Point p2 = pair.second;
    if (point_inside_polygons(p1, -1) || point_inside_polygons(p2, -1)) {
        return true;
    }
    for (const Edge &e : edge) {
        if (segments_intersection(p1, p2, e, p)) {
            return true;
        }
    }
    return false;
}

bool Graph::segment_polygons_intersection(const Point &p1_, const Point &p2_) {
    auto pair = get_close_points(p1_, p2_);
    Point p1 = pair.first;
    Point p2 = pair.second;

    if (point_inside_polygons(p1, -1) || point_inside_polygons(p2, -1)) {
        return true;
    }
    for (const auto &pol : polygons) {
        Point pt1 = pol.Points[0];
        Point pt2;
        Point pt_intersect;
        for (size_t i = 1; i < pol.Points.size(); ++i) {
            pt2 = pol.Points[i];
            if (segments_intersection(p1, p2, pt1, pt2, pt_intersect)) {
                return true;
            }
            pt1 = pt2;
        }
    }
    return false;
}

std::vector<Point> Graph::get_intersection_points(int polygon_id, Point pt_from, Point pt_to) {
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

Point Rectangle::mid_point() const {
    return mid_points(perp1, perp2);
}

Point Rectangle::get_perpendicular_vector() const {
    return get_vector(base1, perp1);
}

Point Rectangle::center() const {
    return mid_points(base1, perp2);
}

double Edge::calculate_cost() const {
    double value = 0.0;
    if (o_type == TRENCH_EDGE) {
        value = getDist() * config.trench_cost;
    } else if(o_type == HDD_EDGE) {
        value = getDist() * config.hdd_cost;
    } else if (o_type == TRENCH_TO_HDD_EDGE) {
        value = config.trench_hdd_cost;
    }
    return value;
}
