#include <vector>
#include <point.h>
#include <cstdio>

class pointset {
    std::vector<point<2>> points;
public:
    using iterator = std::vector<point<2>>::iterator;

    iterator begin();
    iterator end();

    pointset operator&(const pointset& ro) const;
    pointset operator|(const pointset& ro) const;
};