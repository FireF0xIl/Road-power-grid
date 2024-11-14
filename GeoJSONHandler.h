#ifndef HDD_GEOJSONHANDLER_H
#define HDD_GEOJSONHANDLER_H

class Graph;

#include <nlohmann/json.hpp>
#include "Graph.h"

class GeoJSONHandler {
public:
    GeoJSONHandler();
    void load();
    std::vector<nlohmann::json> get_polygons_coordinates() const;
    nlohmann::json get_geojson_reference() const; // temporary

private:
    void fill_features();

    std::vector<nlohmann::json> features;
    nlohmann::json jf;
};

void save_graph(const Graph &graph, const nlohmann::json &reference);


#endif
