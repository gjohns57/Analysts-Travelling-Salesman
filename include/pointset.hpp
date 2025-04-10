#pragma once

#include <vector>
#include "point.hpp"
#include <cstdio>

class pointset {
    std::vector<point<2> > points;
public:
    using iterator = std::vector<point<2> >::iterator;

    iterator begin();
    iterator end();

    void add_point(point<2> p);

    pointset operator&(const pointset& ro) const;
    pointset operator|(const pointset& ro) const;
};
