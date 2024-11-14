#include "GeoJSONHandler.h"
#include "GeoJSONUtility.h"

GeoJSONHandler::GeoJSONHandler() {
    Config::read();
}

void GeoJSONHandler::load() {
    if (Config::input.empty()) {
        return;
    }
    jf = GeoJSONReader(Config::input);
    fill_features();
}

std::vector<nlohmann::json> GeoJSONHandler::get_polygons_coordinates() const {
    std::vector<nlohmann::json> polygons;
    for (auto const &f : features) {
        if (f["geometry"]["type"] == "Polygon") {
            polygons.emplace_back(f["geometry"]["coordinates"]);
        }
    }
    return polygons;
}

void GeoJSONHandler::fill_features() {
    for (const auto &f : jf["features"]) {
        if (f["type"] == "Feature") {
            features.emplace_back(f);
        }
    }
}

nlohmann::json GeoJSONHandler::get_geojson_reference() const {
    return jf;
}

void save_graph(const Graph &graph, const nlohmann::json &reference) {
    GeoJSONWriter(reference, graph.vertex, V_TYPE::TRENCH, "trench_");
    GeoJSONWriter(reference, graph.vertex, V_TYPE::HDD, "hdd_");
    GeoJSONWriter(reference, graph.edge, EDGE_TYPE::TRENCH_EDGE, "trench_");
    GeoJSONWriter(reference, graph.edge, EDGE_TYPE::HDD_EDGE, "hdd_");
    GeoJSONWriter(reference, graph.edge, EDGE_TYPE::TRENCH_TO_HDD_EDGE, "trench_to_hdd_");
}
