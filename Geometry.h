#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <map>
#include <cmath>
//#include "Graph.h"

const double PI = 3.141592653589793238463;
const double EPS = 1E-9;
//const double DENSITY = 30.0;
const double DENSITY = 25.0;

struct Point {
    double x, y; // координаты точки
    Point() : x(0), y(0) {}
    Point(double x, double y) : x(x), y(y) {}

    bool operator<(const Point & p) const {
        return (x < p.x-EPS) || ((std::abs(x-p.x) < EPS) && (y < p.y - EPS));
    }

    bool operator ==(const Point &p) const {
        return (std::abs(x-p.x) < EPS) && (std::abs(y-p.y) < EPS);
    }

    double dist(const Point &other) const {
        return std::sqrt(std::pow(x - other.x, 2) + std::pow(y - other.y, 2));
    }

    void print() const;
};

Point get_vector(const Point &p1, const Point &p2);

Point getPerpendicularVector(const Point &p1, const Point &p2);

Point get_between_points(const Point &p1, const Point &p2, double w1 = 1, double w2 = 1);

Point mid_points(const Point &p1, const Point &p2);

double dist(double x1, double y1, double x2, double y2);

struct Segment {
    Point p, q;
    Segment (Point p, Point q) : p(p), q(q) {}
};

struct Polygon {
    int id;
    std::vector<Point> Points;
    std::vector<Segment> Segments;

    Polygon(): Points() {}
    Polygon(int _id): id(_id), Points() {}

    void add_point(Point &point);
    void add_point(double first, double second);

    void print() const;

    using iterator = std::vector<Point>::iterator;
    using const_iterator = std::vector<Point>::const_iterator;

    iterator begin() noexcept {
        return Points.begin();
    }
    iterator end() noexcept {
        return Points.end();
    }

    const_iterator begin() const noexcept {
        return Points.begin();
    }
    const_iterator end() const noexcept {
        return Points.end();
    }

    std::size_t size() const {
        return Points.size();
    }

    Point &operator[](size_t i) {
        return Points[i];
    }
    Point const &operator[](size_t i) const {
        return Points[i];
    }
};

struct Rectangle {
    Point base1, base2;
    Point perp1, perp2;
    size_t id1, id2;
    Rectangle (const Point &base1, const Point &base2, const Point &perp1, const Point &perp2) :
            base1(base1), base2(base2), perp1(perp1), perp2(perp2), id1(-1), id2(-1)
    {}
    Rectangle (const Point &base1, const Point &base2, const Point &perp1, const Point &perp2, size_t id1, size_t id2) :
            base1(base1), base2(base2), perp1(perp1), perp2(perp2), id1(id1), id2(id2)
    {}

    Point mid_point() const;

    Point center() const;

    Point get_perpendicular_vector() const;

    bool operator <(const Rectangle &other) const {
        return base1 < other.base1 && base2 < other.base2;
    }
    bool operator ==(const Rectangle &other) const {
        return (base1 == other.base1 && base2 == other.base2 && perp1 == other.perp1 && perp2 == other.perp2) ||
        (base1 == other.base2 && base2 == other.base1 && perp1 == other.perp2 && perp2 == other.perp1);
    }
};

bool segments_intersection(const Point &p1, const Point &q1, const Point &p2, const Point &q2,
                           Point &intersection);


bool point_inside_polygon(const Point &point, Polygon polygon);

double find_alpha(const Point &p1, const Point &p2);

double degree_to_radian(double angle);

int segments_splice(const Point &p1, const Point &p2);

std::vector<Point> splice_segment(const Point &p1, const Point &p2);

std::pair<Point, Point> get_close_points(const Point &p1_, const Point &p2_);

std::pair<Point, Point> get_distant_points(const Point &p1_, const Point &p2_);

double get_distance_from_segment(Point target, Point p1, Point p2);

#endif