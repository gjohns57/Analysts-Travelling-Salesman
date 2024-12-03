#include <cstdio>
#include <cstdlib>
#include <climits>
#include <point.h>
#include <vector>
#include <local_flatness.h>

#define NUM_POINTS 4

double rand_float(double min, double max) {
    return (rand() / (double)RAND_MAX) * (max - min) + min;
}

int main(int argc, char** argv) {

    std::vector<point<2>> points;

    point<2> direction;
    direction[0] = rand_float(-10, 10);
    direction[1] = rand_float(-10, 10);

    printf("\\centering\n");
    printf("\\begin{tikzpicture}\n");

    points.push_back(point<2>());
    for (int i = 0; i < 20; i++) {
        point<2> p;
        p[0] = rand_float(-10, 10);
        p[1] = rand_float(-10, 10);
        double min_dist = INT_MAX;

        for (point<2> q : points) {
            double d = (p - q).norm2();
            if (d < min_dist) {
                min_dist = d;
            }
        }

        if (min_dist < 5 || p.norm2() > 10) {
            i--;
            continue;
        }

        printf("\\filldraw [black] (%lf, %lf) circle (5pt);\n", p[0], p[1]);
        points.push_back(p);
    }
    printf("\\filldraw [red] (%lf, %lf) circle (5pt);\n", 0, 0);

    cylinder c = find_thinnest_cylinder(points, point<2>(), 10, 1);
    printf("\\draw [red, thick] (%lf, %lf) -- (%lf, %lf);\n", c.x0[0], c.x0[1], c.x1[0], c.x1[1]);

    printf("\\end{tikzpicture}\n");


    point<2> offset;

    offset[0] = c.width / 2 * (c.x1[1] - c.x0[1]) / c.x1.norm2();
    offset[1] = c.width / 2 * (c.x0[0] - c.x1[0]) / c.x1.norm2();

    printf("\\begin{tikzpicture}\n");

    for (point<2> p : points) {
        printf("\\filldraw [black] (%lf, %lf) circle (5pt);\n", p[0], p[1]);
    }
    printf("\\filldraw [red] (%lf, %lf) circle (5pt);\n", 0, 0);

    printf("\\draw [red, thick] (%lf, %lf) -- (%lf, %lf);\n", c.x0[0], c.x0[1], c.x1[0], c.x1[1]);

    printf("\\draw [red, thick, dotted] (%lf, %lf) -- (%lf, %lf);\n", c.x0[0] + offset[0], c.x0[1] + offset[1], c.x1[0] + offset[0], c.x1[1] + offset[1]);
    printf("\\draw [red, thick, dotted] (%lf, %lf) -- (%lf, %lf);\n", c.x0[0] - offset[0], c.x0[1] - offset[1], c.x1[0] - offset[0], c.x1[1] - offset[1]);
    printf("\\end{tikzpicture}\n");


    return 0;
}