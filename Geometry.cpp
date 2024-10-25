#include <iostream>
#include "Geometry.h"

#include <vector>
#include <cmath>

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
bool segmentsIntersection(const Segment &s1, const Segment &s2,
                          Point &intersection) {
    return segmentsIntersection(s1.p, s1.q, s2.p, s2.q, intersection);
}

bool segmentsIntersection(const Point &p1, const Point &q1, const Point &p2, const Point &q2,
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

double minDistance(const Segment &segment, const Point &E) {
    return minDistance(segment.p, segment.q, E);
}

double radToDegree(double rad) {
    return rad * (180.0/PI);
}

double findAlpha(const Point &p1, const Point &p2) {
    double h1 = hypot(p1.x, p1.y);
    double h2 = hypot(p2.x, p2.y);
    double x1 = p1.x / h1;
    double y1 = p1.y / h1;
    double x2 = p2.x / h2;
    double y2 = p2.y / h2;
    double angle = std::acos(x1 * x2 + y1 * y2);

    if (angle > PI / 2) {
        angle = PI - angle;
    }
    return angle;
}

bool pointInsidePolygon(const Point &point, Polygon polygon) {
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

void Polygon::addPoint(Point &point) {
    this->Points.emplace_back(point);
}

void Polygon::addPoint(double first, double second) {
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

double degreeToRadian(double degree) {
    return (degree * (PI / 180));
}

Point getBetweenPoints(const Point &p1, const Point &p2, double w1, double w2) {
    if (w1 == 0) {
        return p2;
    }
    if (w2 == 0) {
        return p1;
    }
    double s = w1 + w2;
    return {(p1.x * w1 + p2.x * w2) / s, (p1.y * w1 + p2.y * w2) / s};
}

Point midPoints(const Point &p1, const Point &p2) {
    return getBetweenPoints(p1, p2, 1.0, 1.0);
}

Point getVector(const Point &p1, const Point &p2) {
    return {p1.x - p2.x, p1.y - p2.y};
}

Point getPerpendicularVector(const Point &p1, const Point &p2) {
    return {-(p1.y - p2.y), p1.x - p2.x};
}

int segmentsSplice(const Point &p1, const Point &p2) {
    return ceil(p1.dist(p2) / DENSITY);
}


std::vector<Point> spliceSegment(const Point &p1, const Point &p2) {
    std::vector<Point> points{ p1 };
    int s = segmentsSplice(p1, p2);
    for (int i = 1; i <= s; i++) {
        points.emplace_back(getBetweenPoints(p1, p2, s - i, i));
    }
    return points;
}
