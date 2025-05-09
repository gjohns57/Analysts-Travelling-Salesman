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

    for (int i = 0;i<num_points;i++) { // build point<2> vector from Points
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

    unordered_map<pair<Point, Point>, double,PointPairHash,PointPairEqual> distances = generate_distances(Points);
    vector<int> scales = generate_S(Points, distances, v0);

    cout << "Scales:\n";
    for (int scale : scales) {
        cout << scale << '\n';
    }

    set<vector<Point> > nets = generate_all_nets(Points, v0, scales, distances);
    vector<pair<vector<point<2> >,vector<tuple<point<2>,cylinder,vector<point<2> > > > > > flat_points; // first is a net, second is
    // a tuple of a flat point, its cylinder, and the points in its ball
    size_t k = 0; // represents the kth net
    bool val = true;
    vector<vector<point<2> > > p_nets;
    // for loop below builds all the F_k's
    for (vector<Point> net : nets) { // go through each net, we test each point in each net for flatness
        // convert Point net to a point<2> net
        vector<point<2> > p_net;
        for (Point p : net) {
            p_net.push_back(convert_Point(p));
        }
        p_nets.push_back(p_net);
        vector<tuple<point<2>,cylinder,vector<point<2> > > > flat_points_in_net; // each flat point is with its cylinder and ball
        for (point<2> p_2 : p_net) { // to test a point for flatness, we have to run the find_thinnest_cylinder
            // function on the points in the ball of size epsilon
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
        flat_points.push_back(make_pair(p_net,flat_points_in_net));
        if (!val) {
            break;
        }
        // p_net.clear(); don't need to clear p_net since it's shadowed in each loop
    }
    for (pair<vector<point<2> >, vector<tuple<point<2>,cylinder,vector<point<2> > > > > pair_flat_net : flat_points) {
        if (pair_flat_net.first.size() == 1 || pair_flat_net.first.size() == 0) { // error check
            continue;
        }

        vector<pair<point<2>,double> > transformed_vals;
        enum status {
            NONE,
            LEFT,
            RIGHT,
            BOTH
        };
        for (tuple<point<2>,cylinder,vector<point<2> > > flat_point : pair_flat_net.second) {
            for (point<2> p : get<2>(flat_point)) {
                transformed_vals.push_back(make_pair(p,transformation(get<0>(flat_point), get<1>(flat_point),p)));
            }
            // sort the vector so that we can add edges easier
            sort(transformed_vals.begin(),transformed_vals.end(),[](pair<point<2>,double> a, pair<point<2>,double> b) {
                return a.second < b.second; // sorts in ascending order
            });
            // now the transformed_vals vector is filled with pairs of points and their transformed vals
            // todo: traverse the trans_val vect and see what sides the points are on. adding the edges should be easy
            enum status current_status = NONE;
            int flat_point_index = 0;
            for (int i = 0;i<transformed_vals.size();i++) {
                if (transformed_vals[i].second == 0) { // store the flat_point_index in the transformed vals vect
                    flat_point_index = i;
                }
                // we cover all the cases
                if (transformed_vals[i].second < 0 && current_status == NONE) {
                    current_status = LEFT;
                }
                else if (transformed_vals[i].second > 0 && current_status == NONE) {
                    current_status = RIGHT;
                }
                else if (transformed_vals[i].second > 0 && current_status == LEFT) {
                    current_status = BOTH;
                }
                else if (transformed_vals[i].second < 0 && current_status == RIGHT) {
                    current_status = BOTH;
                }
            }
            if (current_status == LEFT) { // go through all the points to the left of v and connect in a chain

                point<2> previous = transformed_vals[0].first; // this is the first point that we'll use to make an edge; it'll also be used to store the previous point
                int previous_index = 0; // we don't actually know this yet, we find it in the following loop
                for (pair<point<2>,int> assignment : assignments) {
                    if (assignment.first[0] == previous[0] && assignment.first[1] == previous[1]) {
                        // this for loop finds the assignment index for the first point in the transformed vals vector
                        previous_index = assignment.second;
                    }
                }
                // starting at 1 since we already found i = 0
                for (int i = 1;i<=flat_point_index;i++) { // walk UP from the left

                    for (pair<point<2>,int> assignment : assignments) { // we now have to find the numerical index of this point in the graph adjmatrix

                        if (assignment.first[0] == get<0>(flat_point)[0] && assignment.first[1] == get<0>(flat_point)[1]) {
                            adjmatrix[assignment.second][previous_index] = 1;
                            adjmatrix[previous_index][assignment.second] = 1;
                            previous_index = assignment.second; // this is the key; we change the last vertex to previous_index to keep the chain going
                        }
                    }
                }
            }

            else if (current_status == RIGHT) {
                // this is the first point that we'll use to make an edge; it'll also be used to store the previous point
                point<2> previous = transformed_vals[transformed_vals.size() - 1].first;
                int previous_index = 0; // we find this in the next for loop
                for (pair<point<2>,int> assignment : assignments) {
                    if (assignment.first[0] == previous[0] && assignment.first[1] == previous[1]) {
                        // this for loop finds the assignment index for the last point in the transformed vals vector
                        previous_index = assignment.second;
                    }
                }
                // starting at transformed_vals.size() - 2 since we already did transformed_vals.size() - 1
                for (int i = transformed_vals.size() - 2;i>=flat_point_index;i--) { // walk DOWN from the right

                    for (pair<point<2>,int> assignment : assignments) { // we now have to find the numerical index of this point in the graph adjmatrix

                        if (assignment.first[0] == get<0>(flat_point)[0] && assignment.first[1] == get<0>(flat_point)[1]) {
                            adjmatrix[assignment.second][previous_index] = 1;
                            adjmatrix[previous_index][assignment.second] = 1;
                            previous_index = assignment.second; // this is the key; we change the last vertex to previous_index to keep the chain going
                        }
                    }
                }
            }
            // BOTH doesn't need anything to be done
        }
    }
    // print nets, flat points, etc.
    int count = 1;
    for (auto entry : flat_points) {
        cout << "Net " << count << ": \n";
        for (point<2> p : entry.first) { // print out nets
            cout << "(" << p[0] << ", " << p[1] << "), ";
        }
        cout << "\nFlat points in net " << count << ":\n";
        if (entry.second.empty()) {
            cout << "No flat points in net.\n";
        }
        for (auto tup : entry.second) {
            if (!get<2>(tup).empty()) {
                cout << "(" << get<0>(tup)[0] << ", " <<get<0>(tup)[1] << ") \n";
            }
        }
        cout << '\n';
        count++;
    }
}
