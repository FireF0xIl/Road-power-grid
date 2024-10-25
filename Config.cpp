#include <fstream>
#include <iostream>
#include "Config.h"
#include "Geometry.h"


void Config::read() {
    std::ifstream ifs(R"(.\..\settings.json)");

    if (ifs.good()) {
        nlohmann::json jfc = nlohmann::json::parse(ifs);

        if (!jfc.contains("input")) {
            return;
        } else {
            input = jfc["input"];
        }

        if (jfc.contains("minHDDDistance") && (jfc["minHDDDistance"].is_number())) {
            min_hdd_distance = jfc["minHDDDistance"];
        } else {
            min_hdd_distance = MIN_HDD_DISTANCE;
        }

        if (jfc.contains("maxHDDDistance") && jfc["maxHDDDistance"].is_number()) {
            max_hdd_distance = jfc["maxHDDDistance"];
        } else {
            max_hdd_distance = MAX_HDD_DISTANCE;
        }

        if (jfc.contains("alpha") && jfc["alpha"].is_number()) {
            alpha = jfc["alpha"];
        } else {
            alpha = ALPHA;
        }

        if (jfc.contains("trenchCost") && jfc["trenchCost"].is_number()) {
            trench_cost = jfc["trenchCost"];
            if (jfc.contains("HDDCost") && jfc["HDDCost"].is_number()) {
                hdd_cost = jfc["HDDCost"];
            } else {
                hdd_cost = trench_cost * 2.0;
            }
        } else {
            if (jfc.contains("HDDCost") && jfc["HDDCost"].is_number()) {
                hdd_cost = jfc["HDDCost"];
                trench_cost = hdd_cost / 2.0;
            } else {
                trench_cost = TRENCH_COST;
                hdd_cost = HDD_COST;
            }
        }

        if (jfc.contains("trenchHDDCost") && jfc["trenchHDDCost"].is_number()) {
            trench_hdd_cost = jfc["trenchHDDCost"];
        } else {
            trench_hdd_cost = TRENCH_HDD_COST;
        }

        if (jfc.contains("maxDistanceFromRoad") && jfc["maxDistanceFromRoad"].is_number()) {
            max_distance_from_road = jfc["maxDistanceFromRoad"];
        } else {
            max_distance_from_road = MAX_DISTANCE_FROM_ROAD;
        }
    } else {
        min_hdd_distance = MIN_HDD_DISTANCE;
        max_hdd_distance = MAX_HDD_DISTANCE;
        alpha = ALPHA;
        trench_cost = TRENCH_COST;
        hdd_cost = HDD_COST;
        trench_hdd_cost = TRENCH_HDD_COST;
        max_distance_from_road = MAX_DISTANCE_FROM_ROAD;
    }

    if (min_hdd_distance < 0.0) {
        min_hdd_distance = 0.0;
    }

    if (max_hdd_distance < 0.0) {
        max_hdd_distance = MAX_HDD_DISTANCE;
    }

    if (trench_cost < 0.0) {
        trench_cost = TRENCH_COST;
    }

    if (hdd_cost < 0.0) {
        hdd_cost = HDD_COST;
    }

    if (trench_hdd_cost < 0.0) {
        trench_hdd_cost = TRENCH_HDD_COST;
    }

    if (max_distance_from_road < 0.0) {
        max_distance_from_road = MAX_DISTANCE_FROM_ROAD;
    }

    if (min_hdd_distance > max_hdd_distance) {
        max_hdd_distance = min_hdd_distance;
    }

    if (alpha < 0) {
        alpha = 0;
    }

    if (alpha > 90) {
        alpha = 90;
    }
    alpha = degreeToRadian(alpha);
    std::cout << "Config: " << min_hdd_distance << " : " << max_hdd_distance << " : " << alpha
        << " : " << trench_cost << " : " << hdd_cost << " : " << trench_hdd_cost << " : " << max_distance_from_road <<'\n';
}
