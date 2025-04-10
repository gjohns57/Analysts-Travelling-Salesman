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

class Point { // poient class that stores X and Y values
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
    size_t operator()(const std::pair<Point, Point>& p) const {
            size_t h1 = std::hash<double>()(p.first.x_value()) ^ (std::hash<double>()(p.first.y_value()) << 1);
            size_t h2 = std::hash<double>()(p.second.x_value()) ^ (std::hash<double>()(p.second.y_value()) << 1);
            return h1 ^ (h2 << 1);
        }
};

struct PointPairEqual { // equality tester for the unordered_map
    bool operator()(const std::pair<Point, Point>& p1, const std::pair<Point, Point>& p2) const {
        return p1.first == p2.first && p1.second == p2.second;
    }
};

std::string output_point(Point& a) { // function for testing and output
    std::stringstream ss;
    ss << "(" << a.x_value() << ", " << a.y_value() << ")";
    return ss.str();
}


std::unordered_map<std::pair<Point, Point>, double,PointPairHash,PointPairEqual> generate_distances(std::vector<Point> points) { // function to create table of distances
    std::unordered_map<std::pair<Point, Point>, double,PointPairHash,PointPairEqual> distances; // initialize table
    for (size_t i = 0;i<points.size();i++) {
        for (size_t j = i;j<points.size();j++) {
            distances.insert(make_pair(std::make_pair(points[i],points[j]),points[i].dist_calc(points[j]))); // insert all
            distances.insert(make_pair(std::make_pair(points[j],points[i]),points[i].dist_calc(points[j]))); // distance is symmetric here
        }
    }
    return distances;
}
// the function below generates S, the set of values that we're going to raise 1/2 to
std::vector<int> generate_S(std::vector<Point>& points, std::unordered_map<std::pair<Point,Point>,double,PointPairHash,PointPairEqual>& distances, Point& v0) {
    // make sure that v0 IS IN vector<Point> points
    std::unordered_map<std::pair<Point,Point>,double,PointPairHash,PointPairEqual>::const_iterator jt; // const bc we only need to read from distances
    std::vector<int> set;
    for (size_t i = 0;i<points.size();i++) {
        if (points[i] == v0) {
            break; // let the program through if v0 is in the vertex set
        }
        if (i == (points.size()-1)) {
            std::cerr << "v0 is not in vertex set." << '\n';
            return set;
        }
    }
    int element;
    for (size_t i = 0;i<points.size();i++) {
        if (v0 == points[i]) {
            continue; // we don't want v0 = v0
        }
        jt = distances.find(std::make_pair(v0,points[i])); // access the key-value pair in distances
        element = floor(-1*(log2((jt->second))));
        for (int j = -2;j<=2;j++) {
            set.push_back(element+j); // add the elements to the new set
        }
    }
    return set;
}

std::set<std::vector<Point> > generate_all_nets(std::vector<Point> points, Point v0, std::vector<int> numbers, std::unordered_map<std::pair<Point,Point>,double,PointPairHash,PointPairEqual> distances) {
    std::set<std::vector<Point> > nets;
    std::vector<Point> net;
    net.push_back(v0);
    nets.insert(net);
    for (int n_k : numbers) {
        std::vector<Point> previous_net = net; // use  last net
        std::vector<Point> new_net = previous_net; // previous net
        for (Point point : points) {
            if (find(previous_net.begin(), previous_net.end(), point) != previous_net.end()) {
                continue; // skip if the point is already in the net
            }

            bool val = true;
            for (Point p_in_net : previous_net) {
                double distance = distances[std::make_pair(p_in_net, point)];
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
