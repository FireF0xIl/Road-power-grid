#include "GeoJSONUtility.h"
#include <fstream>
#include <iostream>
#include <set>

using json = nlohmann::json;

// параметры чтения
nlohmann::json GeoJSONReader(const std::string &pass) {
    std::ifstream ifs(pass);
    json jf = json::parse(ifs);
    return jf;
}
template <typename T>
json featuresJson(const std::set<T> &elements);
template <typename T, typename V>
json featuresJson(const std::set<T> &elements, V o_type);
template <typename T>
json featuresJson(const std::vector<T> &elements);
template <typename T, typename V>
json featuresJson(const std::vector<T> &elements, V o_type);


template <typename T>
void localGeoJSONWriter(const json &js, const std::set<T> &elements, const std::string &name) { //params
    json out;
    out["type"] = js["type"];
    std::string name_ = js["name"];
    if (!name.empty()) {
        name_ += "_" + name;
    } else {
        name_ += "out";
    }
    out["name"] = name_;
    out["crs"] = js["crs"];
    out["features"] = featuresJson(elements);
    std::ofstream output_file(name_ + ".geojson");
    if (!output_file.is_open())  {
        std::cout << "\n Failed to open output file";
    } else {
        output_file << out;
        output_file.close();
    }
}

template <typename T>
void localGeoJSONWriter(const json &js, const std::vector<T> &elements, const std::string &name) { //params
    json out;
    out["type"] = js["type"];
    std::string name_ = js["name"];
    if (!name.empty()) {
        name_ += "_" + name;
    } else {
        name_ += "out";
    }
    out["name"] = name_;
    out["crs"] = js["crs"];
    out["features"] = featuresJson(elements);
    std::ofstream output_file(name_ + ".geojson");
    if (!output_file.is_open())  {
        std::cout << "\n Failed to open output file";
    } else {
        output_file << out;
        output_file.close();
    }
}

template <typename T, typename V>
void localGeoJSONWriter(const json &js, const std::vector<T> &elements, V o_type, const std::string &name) { //params
    json out;
    out["type"] = js["type"];
    std::string name_ = js["name"];
    if (!name.empty()) {
        name_ += "_" + name;
    } else {
        name_ += "out";
    }
    out["name"] = name_;
    out["crs"] = js["crs"];
    out["features"] = featuresJson(elements, o_type);
    std::ofstream output_file(name_ + ".geojson");
    if (!output_file.is_open())  {
        std::cout << "\n Failed to open output file";
    } else {
        output_file << out;
        output_file.close();
    }
}

template <typename T, typename V>
void localGeoJSONWriter(const json &js, const std::set<T> &elements, V o_type, const std::string &name) { //params
    json out;
    out["type"] = js["type"];
    std::string name_ = js["name"];
    if (!name.empty()) {
        name_ += "_" + name;
    } else {
        name_ += "out";
    }
    out["name"] = name_;
    out["crs"] = js["crs"];
    out["features"] = featuresJson(elements, o_type);
    std::ofstream output_file(name_ + ".geojson");
    if (!output_file.is_open())  {
        std::cout << "\n Failed to open output file";
    } else {
        output_file << out;
        output_file.close();
    }
}

void GeoJSONWriter(const json &js, const std::set<Point> &points, const std::string &name) {
    localGeoJSONWriter(js, points, name + "points");
}

void GeoJSONWriter(const json &js, const std::set<Vertex> &vertexes, const std::string &name) {
    localGeoJSONWriter(js, vertexes, name + "vertexes");
}

void GeoJSONWriter(const json &js, const std::vector<Vertex> &vertexes, const std::string &name) {
    localGeoJSONWriter(js, vertexes, name + "vertexes");
}

void GeoJSONWriter(const json &js, const std::vector<Vertex> &vertexes, const V_TYPE o_type, const std::string &name) {
    localGeoJSONWriter(js, vertexes, o_type,  name + "vertexes");
}

void GeoJSONWriter(const json &js, const std::set<Edge> &edges, const std::string &name) {
    localGeoJSONWriter(js, edges, name + "edges");
}

void GeoJSONWriter(const json &js, const std::set<Edge> &edges, const EDGE_TYPE o_type, const std::string &name) {
    localGeoJSONWriter(js, edges, o_type, name + "edges");
}


void GeoJSONWriter(const nlohmann::json &json,
                   const std::set<Point> &points,
                   const std::set<Edge> &edges,
                   const std::string &name) {
    GeoJSONWriter(json, points, name);
    GeoJSONWriter(json, edges, name);

}

json lineJson(const Point &p1, const Point &p2) {
    json js;
    js["type"] = "LineString";
    js["coordinates"] = json::array({
        json::array({p1.x, p1.y}),
        json::array({p2.x, p2.y}),
        });
    return js;
}

json lineJson(const Edge &e) {
    return lineJson(e.point1, e.point2);
}

json pointJson(const Point &point) {
    json js = json::object();
    js["type"] = "Point";
    js["coordinates"] = json::array({point.x, point.y});
    return js;
}

json pointJson(const Vertex &vertex) {
    return pointJson(vertex.position);
}

json geometryFeature(const Point &point) {
    return pointJson(point);
}

json geometryFeature(const Vertex &vertex) {
    return pointJson(vertex);
}

json geometryFeature(const Edge &edge) {
    if (edge.o_type == TRENCH_TO_HDD_EDGE) {
        return pointJson(edge.point1);
    } else {
        return lineJson(edge);
    }
}

template <typename T>
json featuresJson(const std::set<T> &elements) {
    json js = json::array();

    for (const T &e : elements) {
        json cur;
        cur["type"] = "Feature";
        cur["properties"] = json::object();
        cur["geometry"] = geometryFeature(e);
        js.emplace_back(cur);
    }
    return js;
}

template <typename T, typename V>
json featuresJson(const std::set<T> &elements, V o_type) {
    json js = json::array();

    for (const T &e : elements) {
        if (e.o_type == o_type) {
            json cur;
            cur["type"] = "Feature";
            cur["properties"] = json::object();
            cur["geometry"] = geometryFeature(e);
            js.emplace_back(cur);
        }
    }
    return js;
}

template <typename T>
json featuresJson(const std::vector<T> &elements) {
    json js = json::array();

    for (const T &e : elements) {
        json cur;
        cur["type"] = "Feature";
        cur["properties"] = json::object();
        cur["geometry"] = geometryFeature(e);
        js.emplace_back(cur);
    }
    return js;
}

template <typename T, typename V>
json featuresJson(const std::vector<T> &elements, V o_type) {
    json js = json::array();
    for (const T &e : elements) {
        if (e.o_type == o_type) {
            json cur;
            cur["type"] = "Feature";
            cur["properties"] = json::object();
            cur["geometry"] = geometryFeature(e);
            js.emplace_back(cur);
        }
    }
    return js;
}


json lineJson(const Rectangle &r) {
    json js;
    js["type"] = "LineString";
    js["coordinates"] = json::array({
                                     json::array({r.base1.x, r.base1.y}),
                                     json::array({r.base2.x, r.base2.y}),
                                     json::array({r.perp2.x, r.perp2.y}),
                                     json::array({r.perp1.x, r.perp1.y}),
                                     json::array({r.base1.x, r.base1.y})
                                    });
    return js;
}

json featuresJson(const std::vector<Rectangle> &rect) {
    json js = json::array();

    for (const auto &r : rect) {
        json cur;
        cur["type"] = "Feature";
        cur["properties"] = json::object();
        cur["geometry"] = lineJson(r);
        js.emplace_back(cur);
    }
    return js;
}

void GeoJSONWriter(const nlohmann::json &js,
                   const std::vector<Rectangle> &rect,
                   const std::string &name) {
    json out;
    out["type"] = js["type"];
    std::string name_ = js["name"];
    if (!name.empty()) {
        name_ += "_" + name;
    } else {
        name_ += "out";
    }
    out["name"] = name_;
    out["crs"] = js["crs"];
    out["features"] = featuresJson(rect);
    std::ofstream output_file(name_ + "_rectangles.geojson");
    if (!output_file.is_open())  {
        std::cout << "\n Failed to open output file";
    } else {
        output_file << out;
        output_file.close();
    }

}

