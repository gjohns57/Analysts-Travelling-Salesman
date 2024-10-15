#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
using namespace std;

class Graph {
private:
    int V; // number of vertices
    unordered_map<int, vector<pair<int, int> > > adj; // Track edge traversal count
public:
    Graph(int v) : V(v) {}

    void addEdge(int u, int v) {
        adj[u].push_back(make_pair(v, 0));
        adj[v].push_back(make_pair(u, 0));
    }

    vector<int> findEulerTour() {
        vector<int> tour;
        int start_vertex = 1; // Always start at vertex 1

        // Check if vertex 1 exists in the graph
        if (adj.find(1) == adj.end() || adj[1].empty()) {
            cout << "Vertex 1 is not present or is isolated. Cannot start tour from vertex 1." << endl;
            return tour;
        }

        tour.push_back(start_vertex);
        findEulerTourUtil(start_vertex, tour);
        return tour;
    }
    private:
        void findEulerTourUtil(int u, vector<int>& tour) {
            for (vector<pair<int, int> >::iterator it = adj[u].begin(); it != adj[u].end(); ++it) { // iterates through all edges
                int v = it->first; // stores node number
                if (it->second < 2) { // if the node has been traversed less than twice, we go to the next node
                    it->second++; // increment
                    for (vector<pair<int, int> >::iterator jt = adj[v].begin(); jt != adj[v].end(); ++jt) { // iterates through all edges of next vertex
                        if (jt->first == u) { //
                            jt->second++;
                            break;
                        }
                    }
                    findEulerTourUtil(v, tour);
                    tour.push_back(u); // add vertex to
                }
            }
        }
};
int main() {
    int V, E;
    cout << "Enter the number of vertices: ";
    cin >> V;

    Graph graph(V);

    cout << "Enter the number of edges: ";
    cin >> E;

    cout << "Enter the edges (u v):" << endl;
    for (int i = 0; i < E; i++) {
        int u, v;
        cin >> u >> v;
        graph.addEdge(u, v);
    }

    vector<int> tour = graph.findEulerTour();

    if (!tour.empty()) {
        cout << "2-to-1 Euler Tour: ";
        for (vector<int>::const_iterator it = tour.begin(); it != tour.end(); ++it) {
            cout << *it << " ";
        }
        cout << endl;
    }

    return 0;
}
