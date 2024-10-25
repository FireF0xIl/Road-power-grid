#include "main.h"
#include "GeoJSONUtility.h"



int main() {
    Graph g;
    g.load();
    g.build();
    g.save();
    return 0;
}
