#ifndef CONFIG_H
#define CONFIG_H

#include "nlohmann/json.hpp"
#include "Geometry.h"

const constexpr double MIN_HDD_DISTANCE = 30.0;
const constexpr double MAX_HDD_DISTANCE = 150.0;
const constexpr double ALPHA = 10.0;
const constexpr double TRENCH_COST = 1.0;
const constexpr double HDD_COST = 2.0;
const constexpr double TRENCH_HDD_COST = 0.5;
const constexpr double MAX_DISTANCE_FROM_ROAD = 3.0;

class Config {
public:
    static double min_hdd_distance;
    static double max_hdd_distance;
    static double alpha;
    static double trench_cost;
    static double hdd_cost;
    static double trench_hdd_cost;
    static double max_distance_from_road;
    static std::string input;

    static void read();

private:
    static nlohmann::json jf;

};

#endif
