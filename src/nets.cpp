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
};

struct PointPairHash { // hash function for the unordered_map
    size_t operator()(const pair<Point, Point>& p) const {
            size_t h1 = std::hash<double>()(p.first.x_value()) ^ (std::hash<double>()(p.first.y_value()) << 1);
            size_t h2 = std::hash<double>()(p.second.x_value()) ^ (std::hash<double>()(p.second.y_value()) << 1);
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

unordered_map<pair<Point,Point>,double,PointPairHash,PointPairEqual>::const_iterator it; // now I can use this everywhere

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
        element = abs((int)log2((jt->second)));
        for (int j = -2;j<=2;j++) {
            set.push_back(element+j); // add the elements to the new set
        }
    }
    return set;
}

set<vector<Point> > generate_all_nets(vector<Point> points,Point v0,vector<int> numbers,unordered_map<pair<Point,Point>,double,PointPairHash,PointPairEqual> distances) {
    set<vector<Point> > nets; // we use a set of vectors of Points so that we can ignore when the sets come out the same
    vector<Point> net,final_sequence;
    net.push_back(v0); // v0 starts in the net
    for (size_t i = 0;i<numbers.size();i++) { // this loop iterates through set
        for (size_t j = 0;j<points.size();j++) { // this loop generates a net with the power 0.5^set[i]
            vector<Point>::iterator finder = find(net.begin(),net.end(),points[j]);
            if (finder != net.end()) {
                continue; // this means that the element is already in the net
            }
            it = distances.find(make_pair(v0,points[j]));
            if (it->second > pow(0.5,numbers[i])) {
                net.push_back(points[j]);
            }
        }
        nets.insert(net); // store the current net
        // output the net
        cout << "Net " << i << ": ";
        for (size_t output_i = 0;output_i<net.size();output_i++) {
            cout << output_point(net[output_i]) << " ";
        }

        cout << '\n';
    }
    return nets;
}

int main() {
    // initialize everything
    unordered_map<pair<Point, Point>, double, PointPairHash, PointPairEqual> distances;
    double v0x,v0y;
    vector<Point> points;
    set<vector<Point> > nets;
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
    for (size_t i = 0;i<set.size();i++) {
        cout << set[i] << " ";
    }
    cout << '\n';
    nets = generate_all_nets(points, v0, set, distances);
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
