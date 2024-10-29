#include <igraph.h>

#include <point.h>
#include <pointset.h>

#include <vector>
#include <cmath>
#include <iterator>

#define C0 300


double flatness(pointset& ps, point<2> v, double epsilon, double k) {
    int l = 40 * C0;
    double alpha = 1.0e300;;

    for (int i0 = -2 * l; i0 <= 2 * l; i0++) {
        for (int j0 = -2 * l; j0 <= 2 * l; j0++) {
            point<2> x;
            x[0] = C0 * pow(2, -k) * epsilon * (1.0 / l) * i0 + v[0];
            x[1] = C0 * pow(2, -k) * epsilon * (1.0 / l) * j0 + v[1];

            for (int i1 = -2 * l; i1 <= 2 * l; i1++) {
                for (int j1 = -2 * l; j1 <= 2 * l; j1++) {
                    double mindist = 1.0e300;
                    for (point<2> p : ps) {
                        point<2> y;
                        y[0] = C0 * pow(2, -k) * epsilon * (1.0 / l) * i1 + v[0];
                        y[1] = C0 * pow(2, -k) * epsilon * (1.0 / l) * j1 + v[1];
                        if (y[0] != x[0] || y[1] != x[1]) {
                            if (dist(v, line(x, y)) < mindist) {
                                mindist = dist(v, line(x, y));
                            }
                        }
                    }

                    if (mindist < alpha) {
                        alpha = mindist;
                    }
                }
            }
        }
    }

    return alpha;
}