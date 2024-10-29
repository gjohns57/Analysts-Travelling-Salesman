#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
using namespace std;

class Graph {
private:
    int V; // number of vertices
    unordered_map<int, vector<pair<int, int> > > adj; // map to keep track of each vertex, key is vertex number
   // vector of pair of integers is to show where the vertex has edges to, and how many times they've been traversed
public:
    Graph(int v) {
        V = v;
    }
    void addEdge(int u, int v) {
        // add it to both sides;
        adj[u].push_back(make_pair(v, 0));
        adj[v].push_back(make_pair(u, 0));
    }

    vector<int> findEulerTour() {
        vector<int> tour;
        int start_vertex = 1; // always start at vertex 1

        // check if vertex 1 exists in the graph
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
                it->second++; // increment to next node
                for (vector<pair<int, int> >::iterator jt = adj[v].begin(); jt != adj[v].end(); ++jt) { // updates reverse edge
                    if (jt->first == u) { // look for the node that we just came from
                        jt->second++; // increment
                        break; // get out
                    }
                }
                findEulerTourUtil(v, tour); // recursion so we can keep going through the graph until the end
                tour.push_back(u); // add vertex to tour
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
            // we use an iterator here instead of a subscript to iterate through the vector because I didn't make an accessor function
        }
        cout << endl;
    }

    return 0;
}
