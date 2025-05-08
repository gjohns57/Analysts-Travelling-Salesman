#pragma once
#include "point.hpp"
#include <vector>

class cylinder {
public:

    double width;
    point<2> x0, x1;
    cylinder(double width, point<2> x0, point<2> x1) : width(width), x0(x0), x1(x1) {}
};

cylinder find_thinnest_cylinder(std::vector<point<2> > &ps, point<2> v, double epsilon, long k);
