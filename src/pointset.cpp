#include "pointset.h"
#include <algorithm>

void pointset::add_point(point<2> p) {
    points.push_back(p);

    // not super efficient but will have to do
    std::sort(points.begin(), points.end(), [](point<2>& a, point<2>& b) {
        if (a[0] == b[0]) {
            return a[1] < b[1];
        }
        return a[0] < b[0];
    });
}

pointset pointset::operator&(const pointset& ro) const {
    pointset ret;
    std::set_intersection(points.begin(), points.end(), ro.points.begin(), ro.points.end(), std::back_inserter(ret.points));
    return ret;
}

pointset pointset::operator||(const pointset& ro) const {

}