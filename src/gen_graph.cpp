#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <set>
#include <chrono>

using namespace std;

// Structure to represent an edge
struct Edge {
    int v1, v2;

    Edge(int vertex1, int vertex2) {
        // Always store smaller vertex first for consistency
        v1 = min(vertex1, vertex2);
        v2 = max(vertex1, vertex2);
    }

    bool operator<(const Edge& other) const {
        if (v1 != other.v1) return v1 < other.v1;
        return v2 < other.v2;
    }
};

int main() {
    int num_vertices;
    double density;

    // Get input from user
    cout << "Enter number of vertices: ";
    cin >> num_vertices;

    cout << "Enter density (between 0 and 1): ";
    cin >> density;

    if (density < 0 || density > 1) {
        cout << "Density must be between 0 and 1" << endl;
        return 1;
    }

    // Calculate maximum possible edges for an undirected simple graph
    int max_edges = (num_vertices * (num_vertices - 1)) / 2;
    int num_edges = static_cast<int>(density * max_edges);

    // Initialize random number generator
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 gen(seed);
    uniform_int_distribution<> vertex_dist(0, num_vertices - 1);

    // Generate random edges
    set<Edge> edges;
    while (edges.size() < num_edges) {
        int v1 = vertex_dist(gen);
        int v2 = vertex_dist(gen);

        // Skip self-loops
        if (v1 != v2) {
            edges.insert(Edge(v1, v2));
        }
    }

    // Write to file
    ofstream outfile("testing.dim");
    if (!outfile) {
        cout << "Error opening file" << endl;
        return 1;
    }

    // Write number of vertices and edges
    outfile << num_vertices << " " << edges.size() << endl;

    // Write edges
    for (const Edge& edge : edges) {
        outfile << edge.v1 << " " << edge.v2 << endl;
    }

    outfile.close();
    cout << "Graph has been generated and saved to testing.dim" << endl;

    return 0;
}
