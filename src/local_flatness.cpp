#include <local_flatness.h>

#include <vector>
#include <cmath>
#include <iterator>
#include <utility>

#define C0 300

double dist(point<2> p, point<2> a, point<2> b)
{
    // distance from p to the line through a and b
    double x = p[0];
    double y = p[1];
    double x1 = a[0];
    double y1 = a[1];
    double x2 = b[0];
    double y2 = b[1];

    return std::abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / std::sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1));
}

point<2> get_boundary_point(point<2> center, int side, double radius, double offset)
{
    point<2> ret;
    switch (side)
    {
    case 0:
        ret[0] = center[0] + offset;
        ret[1] = center[1] + radius;
        break;
    case 1:
        ret[0] = center[0] + radius;
        ret[1] = center[1] + offset;
        break;
    case 2:
        ret[0] = center[0]+ offset;
        ret[1] = center[1] - radius;
        break;
    case 3:
        ret[0] = center[0] - radius;
        ret[1] = center[1] + offset;
        break;
    }

    return ret;
}

cylinder find_thinnest_cylinder(std::vector<point<2>> &ps, point<2> v, double radius)
{
    int l = 1;
    double width = 1.0e300;
    point<2> cx0, cx1;

    for (int side0 = 0; side0 < 4; side0++)
    {
        for (int side1 = side0 + 1; side1 < 4; side1++)
        {

            for (int i = -2 * l; i <= 2 * l; i++)
            {
                point<2> x0 = get_boundary_point(v, side0, radius, radius * i / (2.0 * l)); // pow(2, -k) * epsilon, C0 * pow(2, -k) * epsilon * (1.0 / l) * i);

                for (int j = -2 * l; j <= 2 * l; j++)
                {
                    point<2> x1 = get_boundary_point(v, side1, radius, radius * j / (2.0 * l)); // pow(2, -k) * epsilon, C0 * pow(2, -k) * epsilon * (1.0 / l) * j);
                    if (x0[0] == x1[0] && x0[1] == x1[1])
                    {
                        continue;
                    }

                    double max_d = 0.0;
                    for (point<2> p : ps)
                    {
                        double d = dist(p, x0, x1);

                        if (d > max_d)
                        {
                            max_d = d;
                        }
                    }

                    if (max_d < width)
                    {
                        width = max_d;
                        cx0 = x0;
                        cx1 = x1;
                    }
                }
            }
        }
    }

    return cylinder(width, cx0, cx1);
}

std::vector<point<2>> get_points_in_ball(point<2> v, std::vector<point<2>> xp, double radius) {
    std::vector<point<2>> ret;

    for (point<2> p : xp) {
        if ((p[0] - v[0]) * (p[0] - v[0]) + (p[1] - v[1]) * (p[1] - v[1]) < radius * radius) {
            ret.push_back(p);
        }
    }

    return ret;
}

std::vector<std::pair<point<2>,point<2>>> flat_pairs(std::vector<point<2>> &x, std::vector<point<2>> &xp, double epsilon, long k) {
    std::vector<std::pair<point<2>,point<2>>> ret;

    for (point<2> p : x) {
        std::vector<point<2>> xp_in_ball = get_points_in_ball(p, xp, epsilon);
        cylinder c = find_thinnest_cylinder(xp_in_ball, p, epsilon);
        double alpha = (1 << k) * c.width / epsilon;

        if (alpha < 1.0 / 16) {
            point<2> next_point;
            double min_cylinar_aligned_component = 1.0e300;
            for (point<2> q : xp_in_ball) {
                if ((q - p).norm2() < epsilon || (q - p).norm2() >= C0 * pow(2, -k - 1) * epsilon) {
                    continue;
                }

                double component = fabs((q - p) * (c.x1 - c.x0) / (c.x1 - c.x0).norm2());  
                if (component < min_cylinar_aligned_component) {
                    min_cylinar_aligned_component = component;
                    next_point = q;
                }
            }

            if ((next_point - p).norm2() > epsilon && (next_point - p).norm2() < 2 * epsilon) {
                ret.push_back(std::make_pair(p, next_point));
            }
        }
    }
}