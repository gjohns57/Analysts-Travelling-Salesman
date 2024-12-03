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
using namespace std;

class Point { // point class that stores X and Y values
    private:
        double x;
        double y;
    public:
        Point(double X, double Y) {
            x = X;
            y = Y;
        }

        double x_value() const {
            return x;
        }
        double y_value() const {
            return y;
        }

        double dist_calc(Point& b) const { // distance from a point to another point
            return sqrt((x - b.x)*(x - b.x) + (y - b.y)*(y - b.y));
        }

        bool operator==(const Point& point) const {
                return x == point.x && y == point.y;
            }
        bool operator!=(const Point& point) const {
                return !(x == point.x && y == point.y);
            }
        bool operator<(const Point& point) const {
                return x < point.x || (x == point.x && y < point.y);
            }
};

struct PointPairHash { // hash function for the unordered_map
    size_t operator()(const pair<Point, Point>& p) const {
            size_t h1 = hash<double>()(p.first.x_value()) ^ (hash<double>()(p.first.y_value()) << 1);
            size_t h2 = hash<double>()(p.second.x_value()) ^ (hash<double>()(p.second.y_value()) << 1);
            return h1 ^ (h2 << 1);
        }
};

struct PointPairEqual { // equality tester for the unordered_map
    bool operator()(const pair<Point, Point>& p1, const pair<Point, Point>& p2) const {
        return p1.first == p2.first && p1.second == p2.second;
    }
};

string output_point(Point& a) { // function for testing and output
    stringstream ss;
    ss << "(" << a.x_value() << ", " << a.y_value() << ")";
    return ss.str();
}


unordered_map<pair<Point, Point>, double,PointPairHash,PointPairEqual> generate_distances(vector<Point> points) { // function to create table of distances
    unordered_map<pair<Point, Point>, double,PointPairHash,PointPairEqual> distances; // initialize table
    for (size_t i = 0;i<points.size();i++) {
        for (size_t j = i;j<points.size();j++) {
            distances.insert(make_pair(make_pair(points[i],points[j]),points[i].dist_calc(points[j]))); // insert all
            distances.insert(make_pair(make_pair(points[j],points[i]),points[i].dist_calc(points[j]))); // distance is symmetric here
        }
    }
    return distances;
}
// the function below generates S, the set of values that we're going to raise 1/2 to
vector<int> generate_S(vector<Point>& points, unordered_map<pair<Point,Point>,double,PointPairHash,PointPairEqual>& distances, Point& v0) {
    // make sure that v0 IS IN vector<Point> points
    unordered_map<pair<Point,Point>,double,PointPairHash,PointPairEqual>::const_iterator jt; // const bc we only need to read from distances
    vector<int> set;
    for (size_t i = 0;i<points.size();i++) {
        if (points[i] == v0) {
            break; // let the program through if v0 is in the vertex set
        }
        if (i == (points.size()-1)) {
            cerr << "v0 is not in vertex set." << endl;
            return set;
        }
    }
    int element;
    for (size_t i = 0;i<points.size();i++) {
        if (v0 == points[i]) {
            continue; // we don't want v0 = v0
        }
        jt = distances.find(make_pair(v0,points[i])); // access the key-value pair in distances
        element = floor(-1*(log2((jt->second))));
        for (int j = -2;j<=2;j++) {
            set.push_back(element+j); // add the elements to the new set
        }
    }
    return set;
}

set<vector<Point> > generate_all_nets(vector<Point> points, Point v0, vector<int> numbers, unordered_map<pair<Point,Point>,double,PointPairHash,PointPairEqual> distances) {
    set<vector<Point> > nets;
    vector<Point> net;
    net.push_back(v0);
    nets.insert(net);
    for (int n_k : numbers) {
        vector<Point> previous_net = net; // use  last net
        vector<Point> new_net = previous_net; // previous net
        for (Point point : points) {
            if (find(previous_net.begin(), previous_net.end(), point) != previous_net.end()) {
                continue; // skip if the point is already in the net
            }

            bool val = true;
            for (Point p_in_net : previous_net) {
                double distance = distances[make_pair(p_in_net, point)];
                if (distance < pow(0.5, n_k)) {
                    val = false;
                    break;
                }
            }

            if (val) {
                new_net.push_back(point);
            }
        }

        net = new_net;
        nets.insert(net);
    }

    return nets;
}


int main() {
    // initialize everything
    unordered_map<pair<Point, Point>, double, PointPairHash, PointPairEqual> distances;
    double v0x,v0y;
    vector<Point> points;
    set<vector<Point> > nets;
    set<vector<Point> >::const_iterator it;
    int num_points;
    double dum_x,dum_y;
    // prompt user
    cout << "How many points?" << endl;
    cin >> num_points;
    cout << "Enter the points: (x y)" << endl;
    for (int i = 0;i<num_points;i++) { // create vertex set
        cin >> dum_x >> dum_y;
        points.push_back(Point(dum_x,dum_y));
    }
    distances = generate_distances(points); // create distances map
    cout << "Enter starting point: (x y)" << endl;
    cin >> v0x >> v0y;
    Point v0 = Point(v0x,v0y);
    vector<int> set = generate_S(points, distances,v0);
    cout << "S = {";
    for (size_t i = 0;i<set.size();i++) {
        cout << set[i] << " ";
    }
    cout << "}" << '\n' << '\n';
    nets = generate_all_nets(points, v0, set, distances);
    int counter = 1;
    for (it = nets.begin();it != nets.end();it++) {
      cout << "Net number " << counter << ": { ";
      for (Point point : *it) {
        cout << output_point(point) << " ";
      }
      cout << "}" << endl;
      counter++;
    }
    // now we have to generate a table of n_k

    /*  output
    for (size_t i = 0;i<points.size();i++) {
        for (size_t j = i;j<points.size();j++) {
            it = distances.find(make_pair(points[i],points[j]));
            string test = output_point(points[i]) + ", " + output_point(points[j]);
            cout << "Distance between " << test << " is " << it->second << endl;
        }
    }
    */
    return 0;
}
