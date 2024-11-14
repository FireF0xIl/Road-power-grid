#include <iostream>
#include "Geometry.h"

#include <vector>
#include <cmath>
#include <numbers>

bool onSegment(Point p, Point q, Point r) {
    if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) && q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y)) {
        return true;
    } else {
        return false;
    }
}

int orientation(Point p, Point q, Point r) {
    double val = (q.y - p.y) * (r.x - q.x)
                 - (q.x - p.x) * (r.y - q.y);
    if (val == 0)
        return 0;
    return (val > 0) ? 1 : 2;
}

bool segments_intersection(const Point &p1, const Point &q1, const Point &p2, const Point &q2,
                           Point &intersection
                          ) {
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);

    if (o1 != o2 && o3 != o4) {
        double a1 = q1.y - p1.y;
        double b1 = p1.x - q1.x;
        double c1 = a1 * p1.x + b1 * p1.y;

        double a2 = q2.y - p2.y;
        double b2 = p2.x - q2.x;
        double c2 = a2 * p2.x + b2 * p2.y;

        double determinant = a1 * b2 - a2 * b1;

        if (determinant != 0) {
            intersection.x
                    = (c1 * b2 - c2 * b1) / determinant;
            intersection.y
                    = (a1 * c2 - a2 * c1) / determinant;
            return true;
        }
    }

    if (o1 == 0 && onSegment(p1, p2, q1))
        return true;
    if (o2 == 0 && onSegment(p1, q2, q1))
        return true;
    return false;
}

double minDistance(const Point &A, const Point &B, const Point &E) {
    double x1 = B.x - A.x;
    double y1 = B.y - A.y;
    double x2 = E.x - A.x;
    double y2 = E.y - A.y;
    return std::abs(x1 * y2 - y1 * x2) / sqrt(x1 * x1 + y1 * y1);
}

double radToDegree(double rad) {
    return rad * (180.0 * std::numbers::inv_pi);
}

double find_alpha(const Point p1, const Point p2) {
    double h1 = hypot(p1.x, p1.y);
    double h2 = hypot(p2.x, p2.y);
    double x1 = p1.x / h1;
    double y1 = p1.y / h1;
    double x2 = p2.x / h2;
    double y2 = p2.y / h2;
    double angle = std::acos(x1 * x2 + y1 * y2);

    if (angle > std::numbers::pi / 2) {
        angle = std::numbers::pi - angle;
    }
    return angle;
}

bool point_inside_polygon(const Point &point, const Polygon &polygon) {
    int num_vertices = polygon.size();
    double x = point.x, y = point.y;
    bool inside = false;

    Point p1 = polygon[0], p2;

    for (int i = 1; i < num_vertices; i++) {
        p2 = polygon[i];

        if (y > std::min(p1.y, p2.y)) {
            if (y <= std::max(p1.y, p2.y)) {
                if (x <= std::max(p1.x, p2.x)) {
                    double x_intersection
                            = (y - p1.y) * (p2.x - p1.x)
                              / (p2.y - p1.y)
                              + p1.x;

                    if (p1.x == p2.x
                        || x <= x_intersection) {
                        inside = !inside;
                    }
                }
            }
        }

        p1 = p2;
    }

    return inside;
}

void Polygon::add_point(Point point) {
    this->Points.emplace_back(point);
}

void Polygon::add_point(double first, double second) {
    this->Points.emplace_back(first, second);
    if (this->Points.size() > 1) {
        Point &p = this->Points[this->Points.size() - 2];
    }
}

void Point::print() const {
    std::cout << "[{" << this->x << "} {}]" << this->y;
}

void Polygon::print() const {
    for (const Point &p : this->Points) {
        p.print();
        std::cout << std::endl;
    }
}

double dist(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x1 - x2, 2) + std::pow(y1 - y2, 2));
}

double degree_to_radian(double degree) {
    return degree * (std::numbers::pi / 180);
}

Point get_between_points(const Point &p1, const Point &p2, double w1, double w2) {
    if (w1 == 0) {
        return p2;
    }
    if (w2 == 0) {
        return p1;
    }
    double s = w1 + w2;
    return {(p1.x * w1 + p2.x * w2) / s, (p1.y * w1 + p2.y * w2) / s};
}

Point mid_points(const Point &p1, const Point &p2) {
    return get_between_points(p1, p2, 1.0, 1.0);
}

Point get_vector(const Point &p1, const Point &p2) {
    return {p1.x - p2.x, p1.y - p2.y};
}

Point getPerpendicularVector(const Point &p1, const Point &p2) {
    return {-(p1.y - p2.y), p1.x - p2.x};
}

int segments_splice(const Point &p1, const Point &p2) {
    return ceil(p1.dist(p2) / DENSITY);
}


std::vector<Point> splice_segment(const Point &p1, const Point &p2) {
    std::vector<Point> points{ p1 };
    int s = segments_splice(p1, p2);
    for (int i = 1; i <= s; i++) {
        points.emplace_back(get_between_points(p1, p2, s - i, i));
    }
    return points;
}

std::pair<Point, Point> get_close_points(const Point &p1_, const Point &p2_) {
    double dx = p2_.x - p1_.x;
    double dy = p2_.y - p1_.y;
    // norm
    double base = hypot(dx, dy);
    dx = dx / base;
    dy = dy / base;
    Point p1 = {p1_.x + dx, p1_.y + dy};
    Point p2 = {p2_.x - dx,  p2_.y - dy};
    return std::make_pair(p1, p2);
}

std::pair<Point, Point> get_distant_points(const Point &p1_, const Point &p2_) {
    double dx = p2_.x - p1_.x;
    double dy = p2_.y - p1_.y;
    // norm
    double base = hypot(dx, dy);
    dx = dx / base;
    dy = dy / base;
    Point p1 = {p1_.x - dx, p1_.y - dy};
    Point p2 = {p2_.x + dx,  p2_.y + dy};
    return std::make_pair(p1, p2);
}


//https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
double get_distance_from_segment(Point target, Point p1, Point p2) {
    double a = target.x - p1.x;
    double b = target.y - p1.y;
    double c = p2.x - p1.x;
    double d = p2.y - p1.y;

    double lenSq = c * c + d * d;

    double param;
    if (lenSq != 0.0) {
        auto dot = a * c + b * d;
        param = dot / lenSq;
    } else {
        param = -1.0;
    }
    double xx;
    double yy;
    if (param < 0) {
        xx = p1.x;
        yy = p1.y;
    } else if (param > 1) {
        xx = p2.x;
        yy = p2.y;
    } else {
        xx = p1.x + param * c;
        yy = p1.y + param * d;
    }
    double dx = target.x - xx;
    double dy = target.y - yy;
    return hypot(dx, dy);
}

Point Rectangle::mid_point() const {
    return mid_points(perp1, perp2);
}

Point Rectangle::get_perpendicular_vector() const {
    return get_vector(base1, perp1);
}

Point Rectangle::center() const {
    return mid_points(base1, perp2);
}
