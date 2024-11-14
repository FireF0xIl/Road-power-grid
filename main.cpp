#include "Graph.h"


int main() {
    Graph g;
    Map map = Map();
    GeoJSONHandler handler = GeoJSONHandler();
    handler.load();
    map.load_polygons(handler);

    g.build(map);
    save_graph(g, handler.get_geojson_reference());
    return 0;
}
