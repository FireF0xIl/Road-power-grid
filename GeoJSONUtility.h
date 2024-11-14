#ifndef GEOJSONUTILITY_H
#define GEOJSONUTILITY_H

#include <nlohmann/json.hpp>
#include <set>
//#include "Geometry.h"
//#include "Graph.h"
#include "GraphElement.h"

nlohmann::json GeoJSONReader(const std::string &pass);


void GeoJSONWriter(const nlohmann::json &json,
                   const std::set<Point> &points,
                   const std::set<Edge> &edges,
                   const std::string &name="");

void GeoJSONWriter(const nlohmann::json &json,
                   const std::set<Point> &points,
                   const std::string &name="");

void GeoJSONWriter(const nlohmann::json &json,
                   const std::set<Vertex> &points,
                   const std::string &name="");

void GeoJSONWriter(const nlohmann::json &json,
                   const std::vector<Vertex> &points,
                   const std::string &name="");

void GeoJSONWriter(const nlohmann::json &json,
                   const std::vector<Vertex> &points,
                   V_TYPE o_type,
                   const std::string &name="");

void GeoJSONWriter(const nlohmann::json &json,
                   const std::set<Edge> &edges,
                   const std::string &name="");

void GeoJSONWriter(const nlohmann::json &json,
                   const std::set<Edge> &edges,
                   EDGE_TYPE o_type,
                   const std::string &name="");

void GeoJSONWriter(const nlohmann::json &json,
                   const std::vector<Rectangle> &rect,
                   const std::string &name="");

#endif