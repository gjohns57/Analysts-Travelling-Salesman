#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <tuple>
#include <utility>
#include "local_flatness_impl.hpp"
#include "nets.hpp"
#include "point.hpp"
using namespace std;

point<2> convert_Point(const Point& p) {
    point<2> result;
    result[0] = p.x_value();
    result[1] = p.y_value();
    return result;
}

Point convert_point(point<2>& p) {
    return Point(p[0],p[1]);
}

point<2> line_proj(point<2>& p, cylinder& c, double& a, double& b) {
    point<2> line_proj;
    line_proj[0] = (p[0] + a*(p[1]-b))/(1+a*a);
    line_proj[1] = (a*((p[0]+a*(p[1]-b))/(1+a*a))+b);
    return line_proj;
}

double transformation(point<2>& v, cylinder& c, point<2>& p) { // cylinder has points for a line
    // calculate slope first
    double a = (c.x0[1]-c.x1[1])/(c.x0[0]-c.x1[0]);
    // we need the arctan of that for the angle in rad
    double angle = 3.1415 - atan(a);
    double b = c.x0[1]-a*c.x1[0];
    // first, we'll project the point onto the line
    point<2> projection = line_proj(p,c,a,b);
    point<2> proj_v = line_proj(v,c,a,b);
    point<2> transformed;
    transformed[0] = (projection[0] - proj_v[0])*cos(angle) - (a*projection[0]+b - proj_v[1])*sin(angle);
    transformed[1] = (projection[0] - proj_v[0])*sin(angle) + (a*projection[0]+b - proj_v[1])*cos(angle);
    return transformed[0];
}

int main() {
    int num_points;
    cout << "How many points would you like to enter? ";
    cin >> num_points;
    double x,y,v0_x,v0_y;
    cout << "Enter your points in the form \"x y\":\n";
    vector<Point> Points;
    vector<point<2> > points;
    vector<vector<int> > adjmatrix(num_points,vector<int>(num_points,0));
    vector<pair<point<2>,int> > assignments;
    for (int i = 0;i<num_points;i++) { // build Points vector
        cin >> x >> y;
        Point p = Point(x,y);
        Points.push_back(p);
    }

    for (int i = 0;i<num_points;i++) { // build points vector from Points
        point<2> temp = convert_Point(Points[i]);
        points.push_back(temp);
        // the assignments are arbitrary, as long as they are consistent
        assignments.push_back(make_pair(temp,i));
    }

    cout << "What are the coordinates of the point you would like to start at? ";
    cin >> v0_x >> v0_y;
    Point v0 = Point(v0_x,v0_y);
    if (find(Points.begin(),Points.end(),v0) == Points.end()) {
        cout << "Starting point was not found in entered points.\n";
        return 1;
    }

    std::unordered_map<std::pair<Point, Point>, double,PointPairHash,PointPairEqual> distances = generate_distances(Points);
    vector<int> scales = generate_S(Points, distances, v0);

    cout << "Scales:\n";
    for (int scale : scales) {
        cout << scale << '\n';
    }

    set<vector<Point> > nets = generate_all_nets(Points, v0, scales, distances);
    vector<pair<vector<point<2> >,vector<tuple<point<2>,cylinder,vector<point<2> > > > > > flat_points; // first is a net, second is
    // a tuple of a flat point, its cylinder, and the points in its ball
    size_t k = 0;
    bool val = true;
    vector<vector<point<2> > > p_nets;
    // for loop below builds all the F_k's
    for (vector<Point> net : nets) { // go through each net, we test each point in each net for flatness
        // convert Point net to a point<2> net
        vector<point<2> > p_net;
        p_nets.push_back(p_net);
        for (Point p : net) {
            p_net.push_back(convert_Point(p));
        }
        vector<tuple<point<2>,cylinder,vector<point<2> > > > flat_points_in_net; // each flat point is with its cylinder and ball
        for (Point p : net) { // to test a point for flatness, we have to run the find_thinnest_cylinder
            // function on the points in the ball of size epsilon
            point<2> p_2 = convert_Point(p);
            if ((k+1) == nets.size()) {
                val = false;
                break;
            }
            double epsilon = 300*pow(0.5,scales[k+1]);
            vector<point<2> > points_in_ball = get_points_in_ball(p_2, p_net, epsilon);
            cylinder c = find_thinnest_cylinder(points_in_ball, p_2, epsilon);
            double alpha = (1 << k) * c.width / epsilon; // should be the flatness
            if (alpha < 1.0/16) {
                flat_points_in_net.push_back(make_tuple(p_2,c,points_in_ball));
            }
        }
        k++;
        if (!flat_points_in_net.empty()) {
            flat_points.push_back(make_pair(p_net,flat_points_in_net));
        }
        if (!val) {
            break;
        }
        p_net.clear(); // clear p_net so we can reuse this vect for the next net
    }
    for (pair<vector<point<2> >, vector<tuple<point<2>,cylinder,vector<point<2> > > > > pair_flat_net : flat_points) {
        if (pair_flat_net.first.size() == 1 || pair_flat_net.first.size() == 0) {
            continue;
        }
        // what do i want to do
        //
        // i need to go through all the flat points in each net, and run the transformation on each point in the net,
        // using the flat point as v. The values from the transformation we'll store as pairs of point<2>s and doubles.
        // Then, we go through that vector and see what side the points are on.
        //
        // if the points are only on one side, then we can go through the other side and connect all the points in the graph
        // if the points are on both sides, then we can stop.
        for (tuple<point<2>,cylinder,vector<point<2> > > flat_point : pair_flat_net.second) {
            vector<pair<point<2>,double> > transformed_vals;
            for (point<2> p : get<2>(flat_point)) {
                transformed_vals.push_back(make_pair(p,transformation(get<0>(flat_point), get<1>(flat_point),p)));
            }
            // now the transformed_vals vector is filled with pairs of points and their transformed vals
        }
    }
    for (pair<vector<point<2> >,vector<tuple<point<2>,cylinder,vector<point<2> >> > > x : flat_points) {
        for (point<2> p : x.first) {
            cout << "Net: \n";
            cout << "(" << p[0] << " " << p[1] << ") ";
        }
        cout << '\n';
        for (tuple<point<2>,cylinder,vector<point<2> > > p : x.second) {
            cout << "Flat points in net: \n";
            cout << "(" << get<0>(p)[0] << " " << get<0>(p)[1] << ") ";
        }
    }
}
