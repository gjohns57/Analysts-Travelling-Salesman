#include <iostream>
#include <iterator>
#include <unordered_set>
#include <utility>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <sstream>
#include <set>
#include <map>
#include "local_flatness.hpp"
#include "nets.hpp"

using namespace std;

class Graph {
private:
    vector<Point> vertices;
    map<Point, vector<Point> > adjacency_list;

public:
    Graph(const vector<Point>& points) {
        vertices = points;
        for (const Point& p : points) {
            adjacency_list[p] = vector<Point>();
        }
    }

    void add_edge(const Point& p1, const Point& p2) {
        adjacency_list[p1].push_back(p2);
        adjacency_list[p2].push_back(p1);
    }

    bool has_edge(const Point& p1, const Point& p2) const {
        if (adjacency_list.find(p1) == adjacency_list.end()) {
            return false;
        }
        const vector<Point>& neighbors = adjacency_list.at(p1);
        return find(neighbors.begin(), neighbors.end(), p2) != neighbors.end();
    }

    void print_graph() const {
        cout << "Graph Representation:" << endl;
        for (const Point& p : vertices) {
            cout << output_point(const_cast<Point&>(p)) << " -> { ";
            const vector<Point>& neighbors = adjacency_list.at(p);
            for (const Point& neighbor : neighbors) {
                cout << output_point(const_cast<Point&>(neighbor)) << " ";
            }
            cout << "}" << endl;
        }
    }
};

point<2> convert_to_point2(const Point& p) {
    point<2> result;
    result[0] = p.x_value();
    result[1] = p.y_value();
    return result;
}

Point convert_to_Point(point<2>& p) { // PROBLEM HERE
    return Point(p[0], p[1]);
}

vector<point<2> > convert_points_vector(const vector<Point>& points) {
    vector<point<2> > result;
    for (const Point& p : points) {
        result.push_back(convert_to_point2(p));
    }
    return result;
}

int main() {
    // Initialize everything
    unordered_map<pair<Point, Point>, double, PointPairHash, PointPairEqual> distances;
    double v0x, v0y;
    vector<Point> points;
    set<vector<Point> > nets;
    int num_points;
    double dum_x, dum_y;

    // Prompt user for input
    cout << "How many points?" << endl;
    cin >> num_points;
    cout << "Enter the points: (x y)" << endl;
    for (int i = 0; i < num_points; i++) {
        cin >> dum_x >> dum_y;
        points.push_back(Point(dum_x, dum_y));
    }

    distances = generate_distances(points);

    cout << "Enter starting point: (x y)" << endl;
    cin >> v0x >> v0y;
    Point v0 = Point(v0x, v0y);

    vector<int> set = generate_S(points, distances, v0);
    cout << "S = {";
    for (size_t i = 0; i < set.size(); i++) {
        cout << set[i] << " ";
    }
    cout << "}" << '\n' << '\n';

    // Generate all nets
    nets = generate_all_nets(points, v0, set, distances);

    // Print the generated nets
    int counter = 1;
    for (const auto& net : nets) {
        cout << "Net number " << counter << ": { ";
        for (Point point : net) {
            cout << output_point(point) << " ";
        }
        cout << "}" << endl;
        counter++;
    }
    Graph graph(points);
    vector<vector<Point> > nets_vector(nets.begin(), nets.end());
    std::set<Point> processed_points;
    for (size_t i = 0; i < nets_vector.size() - 1; i++) {
        vector<Point> current_net = nets_vector[i];
        vector<Point> next_net = nets_vector[i+1];
        // add points from current net to processed
        for (const Point& p : current_net) {
            processed_points.insert(p);
        }
        // find points added in the next net (not in current net)
        vector<Point> new_points;
        for (const Point& p : next_net) {
            if (find(current_net.begin(), current_net.end(), p) == current_net.end()) {
                new_points.push_back(p);
            }
        }
        // convert to point<2> format for local flatness calculations
        vector<point<2> > current_points2 = convert_points_vector(current_net);
        vector<point<2> > next_points2 = convert_points_vector(next_net);

        // calculate flatness between points
        for (size_t j = 0; j < current_net.size(); j++) {
            Point p1 = current_net[j];
            point<2> p1_2 = convert_to_point2(p1);
            // get epsilon value for this net level
            double epsilon = pow(0.5, set[i]);
            long k = set[i];
            // find flat pairs between current point and points in next net
            vector<pair<point<2>, point<2> > > flat_point_pairs = flat_pairs({p1_2}, next_points2, epsilon, k);

            // add edges for flat pairs
            for (const auto& pair : flat_point_pairs) {
                Point flat_p1 = convert_to_Point(pair.first);
                Point flat_p2 = convert_to_Point(pair.second);
                // only add edge if it's not already in the graph
                if (!graph.has_edge(flat_p1, flat_p2)) {
                    graph.add_edge(flat_p1, flat_p2);
                    cout << "Added edge: " << output_point(flat_p1) << " - " << output_point(flat_p2) << endl;
                }
            }
        }
    }
    cout << "Final Graph:" << endl;
    graph.print_graph();
    return 0;
}
