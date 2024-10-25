#ifndef CONFIG_H
#define CONFIG_H

#include "nlohmann/json.hpp"

const double MIN_HDD_DISTANCE = 30.0;
const double MAX_HDD_DISTANCE = 150.0;
const double ALPHA = 10.0;
const double TRENCH_COST = 1.0;
const double HDD_COST = 2.0;
const double TRENCH_HDD_COST = 0.5;
const double MAX_DISTANCE_FROM_ROAD = 3.0;

class Config {
public:
    double min_hdd_distance;
    double max_hdd_distance;
    double alpha;
    double trench_cost;
    double hdd_cost;
    double trench_hdd_cost;
    double max_distance_from_road;
    std::string input;

    void read();

private:
    nlohmann::json jf;

};

#endif
