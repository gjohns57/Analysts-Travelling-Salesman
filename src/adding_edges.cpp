#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <map>
#include <utility>
#include <sstream>
#include <limits>

using namespace std;

class Point {
private:
    double x;
    double y;
public:
    Point(double X = 0.0, double Y = 0.0) : x(X), y(Y) {}

    double x_value() const {
        return x;
    }
    double y_value() const {
        return y;
    }
    double dist_calc(const Point& b) const {
        double dx = x - b.x;
        double dy = y - b.y;
        return sqrt(dx * dx + dy * dy);
    }
    double norm2() const {
        return x * x + y * y;
    }
    double norm() const {
      return sqrt(norm2());
    }

    // operators
    bool operator==(const Point& point) const {
        double epsilon = 1e-9;
        return std::abs(x - point.x) < epsilon && std::abs(y - point.y) < epsilon;
    }
    bool operator!=(const Point& point) const {
        return !(*this == point);
    }
    bool operator<(const Point& point) const {
        double epsilon = 1e-9;
        if (std::abs(x - point.x) > epsilon) {
             return x < point.x;
        }
        return y < point.y - epsilon;
    }

    Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y);
    }
    // Dot product (using Point as a vector)
    double dot(const Point& other) const {
        return x * other.x + y * other.y;
    }
};
struct PointHash {
    size_t operator()(const Point& p) const {
        size_t h1 = hash<double>()(p.x_value());
        size_t h2 = hash<double>()(p.y_value());
        return h1 ^ (h2 << 1);
    }
};
struct PointEqual {
     bool operator()(const Point& p1, const Point& p2) const {
        return p1 == p2;
     }
};

struct PointPairHash {
    size_t operator()(const pair<Point, Point>& p) const {
        PointHash hash_point;
        size_t h1 = hash_point(p.first);
        size_t h2 = hash_point(p.second);
        return h1 ^ (h2 << 1); // Combine hashes
    }
};
struct PointPairEqual {
    bool operator()(const pair<Point, Point>& p1, const pair<Point, Point>& p2) const {
        return p1.first == p2.first && p1.second == p2.second; // just in case
    }
};
string output_point(const Point& a) {
    stringstream ss;
    ss << "(" << a.x_value() << ", " << a.y_value() << ")";
    return ss.str();
}
unordered_map<pair<Point, Point>, double, PointPairHash, PointPairEqual> generate_distances(const vector<Point>& points) {
    unordered_map<pair<Point, Point>, double, PointPairHash, PointPairEqual> distances;
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i; j < points.size(); ++j) {
            double dist = points[i].dist_calc(points[j]);
            distances[make_pair(points[i], points[j])] = dist;
            distances[make_pair(points[j], points[i])] = dist;
        }
    }
    return distances;
}

vector<int> generate_S(const vector<Point>& points, const unordered_map<pair<Point, Point>, double, PointPairHash, PointPairEqual>& distances, const Point& v0) {
    vector<int> s_values;
    bool v0_found = false;
    for(const auto& p : points) {
        if (p == v0) {
            v0_found = true;
            break;
        }
    }
    if (!v0_found) {
         cerr << "Error: Starting point v0 is not in the provided set of points." << endl;
         return s_values; // Return empty vector on error
    }

    for (const auto& point : points) {
        if (v0 == point) {
            continue; // Skip distance to itself
        }
        auto it = distances.find(make_pair(v0, point));
        if (it != distances.end() && it->second > 1e-9) {
             double dist = max(it->second, 1e-12);
            int element = floor(-1.0 * log2(dist));
            for (int j = -2; j <= 2; ++j) {
                s_values.push_back(element + j);
            }
        } else if (it == distances.end()) {
             cerr << "Warning: Distance not found for pair " << output_point(v0) << ", " << output_point(point) << endl;
        }
    }
    sort(s_values.begin(), s_values.end());
    s_values.erase(unique(s_values.begin(), s_values.end()), s_values.end());

    return s_values;
}

vector<vector<Point>> generate_all_nets(const vector<Point>& points, const Point& v0, const vector<int>& numbers, const unordered_map<pair<Point, Point>, double, PointPairHash, PointPairEqual>& distances) {
    vector<vector<Point>> net_sequence;
    vector<Point> current_net;
    bool v0_found = false;
    for(const auto& p : points) {
        if (p == v0) {
            v0_found = true;
            break;
        }
    }
    if (!v0_found) {
         cerr << "Error: Starting point v0 not found in points for net generation." << endl;
         return net_sequence; // ret empty sequence
    }


    current_net.push_back(v0);
    net_sequence.push_back(current_net); // add the initial net {v0}

    for (int n_k : numbers) {
        vector<Point> next_net = current_net; // start with points from the previous net
        double threshold = pow(0.5, n_k);

        for (const Point& point : points) {
            bool already_in_net = false;
            for(const Point& p_in_net : current_net) {
                if (point == p_in_net) {
                    already_in_net = true;
                    break;
                }
            }
            if (already_in_net) {
                continue;
            }

            bool far_enough = true;
            for (const Point& p_in_net : current_net) {
                auto it = distances.find(make_pair(p_in_net, point));
                if (it != distances.end()) {
                    if (it->second < threshold) {
                        far_enough = false;
                        break;
                    }
                } else {
                     cerr << "Error: Distance lookup failed in generate_all_nets." << endl;
                     far_enough = false;
                     break;
                }
            }

            if (far_enough) {
                next_net.push_back(point);
            }
        }
        sort(next_net.begin(), next_net.end()); // sort for consistency, shouldn't affect performance much
        if (next_net.size() > current_net.size()) {
             bool changed = next_net.size() != current_net.size();
              if(!changed){
                for(size_t i = 0; i < next_net.size(); ++i) {
                    if(next_net[i] != current_net[i]) {
                        changed = true;
                        break;
                    }
                }
              }

            if(changed){
                current_net = next_net;
                net_sequence.push_back(current_net);
             }

        }
    }
    return net_sequence;
}

const double C0 = 300.0;

double dist_point_line(const Point& p, const Point& a, const Point& b) {
    double x = p.x_value();
    double y = p.y_value();
    double x1 = a.x_value();
    double y1 = a.y_value();
    double x2 = b.x_value();
    double y2 = b.y_value();

    double normalLength = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
    if (normalLength < 1e-9) {
        return p.dist_calc(a);
    }
    return std::abs((y2 - y1) * x - (x2 - x1) * y + x2 * y1 - y2 * x1) / normalLength;
}

Point get_boundary_point(const Point& center, int side, double radius, double offset) {
    Point ret;
    double cx = center.x_value();
    double cy = center.y_value();
    switch (side) {
    case 0: // Top
        ret = Point(cx + offset, cy + radius);
        break;
    case 1: // Right
        ret = Point(cx + radius, cy + offset);
        break;
    case 2: // Bottom
        ret = Point(cx + offset, cy - radius);
        break;
    case 3: // Left
       ret = Point(cx - radius, cy + offset);
        break;
    default:
        // Should not happen
        ret = center;
        break;
    }
    return ret;
}

class Cylinder {
public:
    double width;
    Point x0, x1;
    Cylinder(double w, Point p0, Point p1) : width(w), x0(p0), x1(p1) {}
};

Cylinder find_thinnest_cylinder(const vector<Point>& ps, const Point& v, double radius) {
    int l = 1;
    double min_width = numeric_limits<double>::max();;
    Point best_x0, best_x1;

    if (ps.empty()) {
       return Cylinder(0.0, v, v); // zero-width cylinder at v
    }


    for (int side0 = 0; side0 < 4; side0++) {
        for (int side1 = side0 + 1; side1 < 4; side1++) {
            for (int i = -2 * l; i <= 2 * l; i++) {
                Point x0 = get_boundary_point(v, side0, radius, radius * i / (2.0 * l));
                for (int j = -2 * l; j <= 2 * l; j++) {
                    Point x1 = get_boundary_point(v, side1, radius, radius * j / (2.0 * l));

                    if (x0 == x1) {
                        continue;
                    }

                    double max_d = 0.0;
                    for (const Point& p : ps) {
                        double d = dist_point_line(p, x0, x1);
                        if (d > max_d) {
                            max_d = d;
                        }
                    }

                    if (max_d < min_width) {
                        min_width = max_d;
                        best_x0 = x0;
                        best_x1 = x1;
                    }
                }
            }
        }
    }
     if (min_width == numeric_limits<double>::max()) {
          return Cylinder(0.0, v, v);
     }
    return Cylinder(min_width, best_x0, best_x1);
}

vector<Point> get_points_in_ball(const Point& v, const vector<Point>& xp, double radius) {
    vector<Point> ret;
    double radius_sq = radius * radius;
    for (const Point& p : xp) {
        double dx = p.x_value() - v.x_value();
        double dy = p.y_value() - v.y_value();
        if (dx * dx + dy * dy < radius_sq) {
            ret.push_back(p);
        }
    }
    return ret;
}

vector<pair<Point, Point>> flat_pairs(const vector<Point>& x, const vector<Point>& xp, double epsilon, long k_param) {
    vector<pair<Point, Point>> ret;
    double epsilon_sq = epsilon * epsilon;

    for (const Point& p : x) {
        vector<Point> xp_in_ball = get_points_in_ball(p, xp, epsilon);

        if (xp_in_ball.empty()) {
            continue;
        }
        Cylinder c = find_thinnest_cylinder(xp_in_ball, p, epsilon);

        double alpha = pow(2.0, k_param) * c.width / epsilon;

        if (alpha < (1.0 / 16.0)) {
            Point next_point = p; // init to p, indicating not found yet
            double min_cylinar_aligned_component = numeric_limits<double>::max();
             double lower_bound_dist_sq = C0 * pow(0.5, k_param + 1) * epsilon; // C0 * 2^-(k+1) * epsilon
             lower_bound_dist_sq *= lower_bound_dist_sq; // Square it

             // vect representing cylinder axis direction
             Point axis_vec = c.x1 - c.x0;
             double axis_vec_norm_sq = axis_vec.norm2();


            for (const Point& q : xp_in_ball) {
                 Point q_minus_p = q - p;
                 double dist_pq_sq = q_minus_p.norm2();
                 if (dist_pq_sq < lower_bound_dist_sq) {
                     continue;
                 }
                 double component = 0.0;
                 if (axis_vec_norm_sq > 1e-12) {
                    // The original code calculates: fabs((q - p) * (c.x1 - c.x0) / (c.x1 - c.x0).norm2());
                    // This looks like the projection length calculation *if* axis_vec was normalized, but uses norm2.
                    // This should calculate projection length correctly: |(q-p) . (axis_vec)| / ||axis_vec||
                     component = std::abs(q_minus_p.dot(axis_vec)) / sqrt(axis_vec_norm_sq);

                 } // else component remains 0.0


                if (component < min_cylinar_aligned_component) {
                    min_cylinar_aligned_component = component;
                    next_point = q;
                }
            }
            if (next_point != p) { // Check if a suitable point q was actually found
                 ret.push_back(make_pair(p, next_point));
            }
        }
    }
    return ret;
}

int main() {
    int num_points;
    cout << "How many points?" << endl;
    cin >> num_points;

    if (num_points <= 0) {
        cerr << "Please enter a positive number of points." << endl;
        return 1;
    }

    vector<Point> points(num_points);
    cout << "Enter the points (x y format, one per line):" << endl;
    for (int i = 0; i < num_points; ++i) {
        double x, y;
        cin >> x >> y;
        points[i] = Point(x, y);
    }

    double v0x, v0y;
    cout << "Enter starting point (x y format):" << endl;
    cin >> v0x >> v0y;
    Point v0(v0x, v0y);

    cout << "\nGenerating distances..." << endl;
    unordered_map<pair<Point, Point>, double, PointPairHash, PointPairEqual> distances = generate_distances(points);

    cout << "Generating S set..." << endl;
    vector<int> s_values = generate_S(points, distances, v0);
    if (s_values.empty()) {
         bool v0_found = false;
         for(const auto& p : points) if(p==v0) v0_found = true;
         if (!v0_found) cerr << "Starting point v0 was not in the list of points." << endl;
         else cerr << "Could not generate S values (maybe v0 is isolated?)." << endl;
         return 1;
    }

    cout << "S = { ";
    for (size_t i = 0; i < s_values.size(); ++i) {
        cout << s_values[i] << (i == s_values.size() - 1 ? "" : ", ");
    }
    cout << " }" << endl;

    cout << "Generating nets..." << endl;
    vector<vector<Point>> net_sequence = generate_all_nets(points, v0, s_values, distances);

    if (net_sequence.empty()) {
        cerr << "Failed to generate any nets." << endl;
        return 1;
    }

    cout << "\nGenerated Net Sequence:" << endl;
    for (size_t i = 0; i < net_sequence.size(); ++i) {
        cout << "Net " << i << " (Size " << net_sequence[i].size() << "): { ";
        for (size_t j = 0; j < net_sequence[i].size(); ++j) {
            cout << output_point(net_sequence[i][j]) << (j == net_sequence[i].size() - 1 ? "" : ", ");
        }
        cout << " }" << endl;
    }

    cout << "\nBuilding graph based on local flatness..." << endl;
    map<Point, vector<Point>> adj_list;
    set<pair<Point, Point>> edges;
    for (size_t k = 0; k < net_sequence.size() - 1; ++k) {
         if (k >= s_values.size()) {
            cerr << "Warning: More nets generated than s_values available. Stopping graph construction." << endl;
             break;
         }
        const vector<Point>& current_net = net_sequence[k];
        const vector<Point>& next_net = net_sequence[k + 1];
        int s_k = s_values[k]; // The exponent used for this transition
        double epsilon_k = pow(0.5, s_k);

        cout << "Analyzing transition Net " << k << " -> Net " << (k + 1) << " with s_k = " << s_k << ", epsilon = " << epsilon_k << endl;
        vector<pair<Point, Point>> pairs_at_step = flat_pairs(current_net, next_net, epsilon_k, s_k);

        cout << " Found " << pairs_at_step.size() << " flat pairs at this step." << endl;

        // Add these pairs as edges to the graph
        for (const auto& pair : pairs_at_step) {
            Point p1 = pair.first;
            Point p2 = pair.second;

            if (p2 < p1) swap(p1, p2);

            if (edges.find(make_pair(p1, p2)) == edges.end()) {
                edges.insert(make_pair(p1, p2));
                adj_list[p1].push_back(p2);
                adj_list[p2].push_back(p1);
                 cout << "  Adding edge: " << output_point(p1) << " -- " << output_point(p2) << endl;
            }
        }
    }
    cout << "\n--- Final Graph Edges ---" << endl;
    if (edges.empty()) {
        cout << "No edges were added to the graph based on the flatness criteria." << endl;
    } else {
        for (const auto& edge : edges) {
            cout << output_point(edge.first) << " -- " << output_point(edge.second) << endl;
        }
    }
     cout << "-------------------------" << endl;
    return 0;
}
