#include <iostream>
#include "Graph.h"
#include <stack>
#include <utility>

void Graph::build(const Map &map) {
    build_trench(map);
    build_hdd(map);
}

void Graph::build_trench(const Map &map) {
    std::cout << "Trench build start" << std::endl;
    for (const auto& ps: map.polygons) {
        Point pt_from = ps.Points[ps.Points.size() - 1];

        if (!map.point_inside_polygons(pt_from, ps.id)) {
            auto result = point.insert(pt_from);
            if (result.second) {
                point_polygon_map.emplace(pt_from, std::set<int>({ps.id}));
            } else {
                point_polygon_map[pt_from].insert(ps.id);
            }
        }

        for (size_t i = 0; i < ps.Points.size(); ++i) {
            Point pt_to = ps.Points[i];
            std::vector<Point> pt_intersect = map.get_intersection_points(ps.id, pt_from, pt_to);
            std::sort(pt_intersect.begin(), pt_intersect.end());
            if (!pt_intersect.empty()) {
                if (ps.Points[i] < pt_from) {
                    std::reverse(pt_intersect.begin(), pt_intersect.end());
                }
                for (const auto pt_i: pt_intersect) {
                    add_trench_point(map, ps.id, pt_from, pt_i);
                    pt_from = pt_i;
                }
                if (!map.point_inside_polygons(pt_to, ps.id)) {
                    add_trench_point(map, ps.id, pt_from, pt_to);
                } else {
                }
            } else {
                if (!map.point_inside_polygons(pt_to, ps.id)) {
                    add_trench_point(map, ps.id, pt_from, pt_to);
                }
            }
            pt_from = ps.Points[i];
        }
    }
    for (auto e: edge) {
        add_edge(V_TYPE::TRENCH, EDGE_TYPE::TRENCH_EDGE, e);
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
    }
    return v_id;
}

void Graph::add_trench_point(const Map &map, int polygon_id, Point pt_from, Point pt_to) {
    if (map.point_inside_polygons(mid_points(pt_from, pt_to), polygon_id)) {
        auto result = point.insert(pt_to);
        if (result.second) {
            point_polygon_map.emplace(pt_to, std::set<int>{polygon_id});
        } else {
            point_polygon_map[pt_to].insert(polygon_id);
        }
    } else {
        auto pt = splice_segment(pt_from, pt_to);
        Point pt_cur = pt_from;
        for (size_t i = 1; i < pt.size(); ++i) {
            auto result = point.insert(pt[i]);
            if (result.second) {
                point_polygon_map.emplace(pt[i], std::set<int>{polygon_id});
            } else {
                point_polygon_map[pt[i]].insert(polygon_id);
            }
            edge.emplace(Edge(EDGE_TYPE::TRENCH_EDGE, pt_cur, pt[i]));
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
    dx = dx / base * 0.5;
    dy = dy / base * 0.5;
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

int choose_sub_graph(const FakeGraph &graph, std::vector<int> &visited) {
    std::map<int, int> col;
    int colour = 0;
    int max = std::numeric_limits<int>::min();
    int chosen_colour = -1;
    for (size_t i = 0; i < graph.mp.size(); i++) {
        if (visited[i] == -1) {
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
                    if (visited[j] == -1) {
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
    return chosen_colour;
}

void Graph::check_hdd_connections(const FakeGraph &graph) {
    std::vector<int> visited(graph.mp.size(), -1);
    int chosen_colour = choose_sub_graph(graph, visited);
    std::vector<size_t> hdd;
    for (size_t i = 0; i < visited.size(); i++) {
        if (visited[i] == chosen_colour) {
            hdd.emplace_back(i);
        }
    }

    std::map<size_t, size_t> translation; // fake to real
    for (auto i : hdd) {
        translation[graph.mp[i].id] = get_vertex_id(V_TYPE::HDD, graph.mp[i].point);
    }
    for (size_t i : hdd) {
        size_t id = translation[i];
        auto cur_graph_mp = graph.mp[i];
        for (auto e : cur_graph_mp.next) {
            double dist = graph.mp[e].point.dist(cur_graph_mp.point);
            double cost = dist * Config::hdd_cost;
            EDGE_TYPE type = EDGE_TYPE::HDD_EDGE;
            vertex[id].graph_edges.emplace_back(type, translation[e], cost);
        }
        for (auto e : cur_graph_mp.real_next) {
            EDGE_TYPE type = EDGE_TYPE::TRENCH_TO_HDD_EDGE;
            double cost = Config::trench_hdd_cost;
            vertex[e].graph_edges.emplace_back(type, id, cost);
            vertex[id].graph_edges.emplace_back(type, e, cost);
        }
    }

    for (size_t i = 0; i < vertex.size(); i++) {
        auto q = vertex[i];
        if (vertex[i].o_type == V_TYPE::HDD) {
            for (auto &e : vertex[i].graph_edges) {
                if (e.vertex_id > i && e.o_type == EDGE_TYPE::HDD_EDGE) {
                    edge.emplace(Edge(EDGE_TYPE::HDD_EDGE, vertex[i].position, vertex[e.vertex_id].position));
                }
            }
        }
    }

    for (size_t i = 0; i < vertex.size(); i++) {
        auto q = vertex[i];
        if (vertex[i].o_type == V_TYPE::TRENCH) {
            for (auto &e : vertex[i].graph_edges) {
                if (e.vertex_id > i && e.o_type == EDGE_TYPE::TRENCH_EDGE) {
                    edge.emplace(Edge(EDGE_TYPE::TRENCH_EDGE, vertex[i].position, vertex[e.vertex_id].position));
                }
            }
        }
    }
    for (size_t i = 0; i < vertex.size(); i++) {
        auto q = vertex[i];
        if (vertex[i].o_type == V_TYPE::TRENCH || vertex[i].o_type == V_TYPE::HDD) {
            for (auto &e : vertex[i].graph_edges) {
                if (e.vertex_id > i && e.o_type == EDGE_TYPE::TRENCH_TO_HDD_EDGE) {
                    edge.emplace(Edge(EDGE_TYPE::TRENCH_TO_HDD_EDGE, vertex[i].position, vertex[e.vertex_id].position));
                }
            }
        }
    }
}

void write_vertex(FakeGraph &local_graph, Point &point, int k, size_t h_id, size_t road1, size_t road2) {
    if (k == 0) {
        local_graph.mp[h_id].real_next.insert(road1);
    } else {
        local_graph.mp[h_id].real_next.insert(road2);
    }
}

void add_hdd_edge(FakeGraph &local_graph,
                  Point &p1,
                  Point &p2,
                  const Rectangle &r1,
                  const Rectangle &r2,
                  int k2,
                  int k4) {
    auto h1_id = local_graph.get_vertex(p1);
    auto h2_id = local_graph.get_vertex(p2);
    local_graph.add_edge(h1_id, h2_id);

    size_t road1 = r1.id1;
    size_t road2 = r1.id2;
    size_t road3 = r2.id1;
    size_t road4 = r2.id2;
    write_vertex(local_graph, p1, k2, h1_id, road1, road2);
    write_vertex(local_graph, p2, k4, h2_id, road3, road4);
}

void Graph::build_hdd(const Map &map) {
    std::cout << "HDD build start" << std::endl;
    std::vector<Rectangle> rectangles;
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
                if (!map.point_inside_polygons(pair.first.mid_point(), -1)) {
                    rectangles.emplace_back(pair.first);
                }
                if (!map.point_inside_polygons(pair.second.mid_point(), -1)) {
                    rectangles.emplace_back(pair.second);
                }
            }
        }
    }

    FakeGraph local_graph = FakeGraph();
    for (int i = 0; i < rectangles.size(); i++) {
        auto r1 = rectangles[i];
        for (int j = i + 1; j < rectangles.size(); j++) {
            auto r2 = rectangles[j];
            double rect_dist = r1.center().dist(r2.center());
            if (rect_dist <= Config::min_hdd_distance - DENSITY || Config::max_hdd_distance + DENSITY <= rect_dist) {
                continue;
            }
            Point pp1 = r1.base1;
            Point pp2 = r1.base2;
            for (int k1 = 0; k1 <= 1; k1++) {
                Point pp3 = r2.base1;
                Point pp4 = r2.base2;
                for (int k2 = 0; k2 <= 1; k2++) {
                    Point p1 = get_between_points(pp1, pp2, 1 - k1, k1);
                    Point perp1;
                    if (k1 == 0) {
                        perp1 = get_between_points(r1.perp1, r1.perp2, 5, 1);
                    } else {
                        perp1 = get_between_points(r1.perp1, r1.perp2, 1, 5);
                    }
                    Point perp2;
                    if (k2 == 0) {
                        perp2 = get_between_points(r2.perp1, r2.perp2, 5, 1);
                    } else {
                        perp2 = get_between_points(r2.perp1, r2.perp2, 1, 5);
                    }
                    Point p2 = get_between_points(pp3, pp4, 1 - k2, k2);
                    if (p1 == p2) {
                        continue;
                    }
                    double dist = p1.dist(p2);
                    if (Config::min_hdd_distance <= dist && dist <= Config::max_hdd_distance) {
                        Point p1p2 = get_vector(p1, p2);
                        auto vector = map.find_segment_intersects(perp1, perp2);
                        if (vector.empty()) {
                            if (map.check_edge_distance(p1, p2)) {
                                add_hdd_edge(local_graph, p1, p2, r1, r2, k1, k2);
                            }
                        } else {
                            if (vector.size() == 1) {
                                continue;
                            }
                            double min_dist = std::numeric_limits<double>::max();
                            double max_dist = std::numeric_limits<double>::min();
                            for (const auto &item : vector) {
                                double cur_dist = item.second.dist(p1);
                                if (cur_dist < min_dist) {
                                    min_dist = cur_dist;
                                }
                                if (cur_dist > max_dist) {
                                    max_dist = cur_dist;
                                }
                            }
                            bool p1_angle = false;
                            bool p2_angle = false;
                            for (const auto &item : vector) {
                                double cur_dist = item.second.dist(p1);
                                if (fabs(cur_dist - min_dist) < EPS) {
                                    if (map.get_alpha(p1p2, item.first) < Config::alpha) {
                                        p1_angle = true;
                                    }
                                }
                                if (fabs(cur_dist - max_dist) < EPS) {
                                    if (map.get_alpha(p1p2, item.first) < Config::alpha) {
                                        p2_angle = true;
                                    }
                                }
                            }
                            if (p1_angle && p2_angle) {
                                if (map.check_edge_distance(p1, p2)) {
                                    add_hdd_edge(local_graph, p1, p2, r1, r2, k1, k2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    check_hdd_connections(local_graph);
    std::cout << "HDD build finish" << std::endl;
}
