#include <iostream>
#include <vector>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <string>
#include <sstream>

typedef std::pair<std::vector<int>, std::vector<std::pair<int, int>>> graph;

void printComponents(std::vector<graph> G);

std::vector<graph> findComponents(std::vector<int> V, std::vector<std::pair<int, int>> E, int node, int *index) {
    std::vector<graph> G;
    std::vector<int> verts;
    std::vector<std::pair<int, int>> edges;
    std::queue<int> q;
    std::unordered_set<int> visited; 
    visited.reserve(V.size());

    int graph_index = 0;

    for(size_t v = 0; v < V.size(); v++) {
        if(visited.find(V[v]) != visited.end()) {
            continue;
        }
        
        q.push(V[v]);
        
        while(!q.empty()) {
            int current = q.front();
            q.pop();

            if(visited.find(current) != visited.end()) {
                continue;
            }
            verts.push_back(current);

            if(current == node) {
                *index = graph_index;
            }

            for(size_t i = 0; i < E.size(); i++) {
                if(current == E[i].first) {
                    q.push(E[i].second);
                    if(std::find(edges.begin(), edges.end(), E[i]) == edges.end()) {
                        edges.push_back(E[i]);
                    }
                }
                else if(current == E[i].second) {
                    q.push(E[i].second);
                    if(std::find(edges.begin(), edges.end(), E[i]) == edges.end()) {
                        edges.push_back(E[i]);
                    }
                }
            }
            visited.insert(current);
        }
        graph g = std::make_pair(verts, edges);
        G.push_back(g);
        graph_index++;

        verts.clear();
        edges.clear();
    }

    printComponents(G);

    return G;
}



graph connectGraph(std::vector<graph> G, int node, int index) {
    size_t total_edge_size = 0;
    size_t total_vert_size = 1;
    if(index == -1) {
        for(size_t i = 0; i < G.size(); i++) {
            G[i].second.push_back(std::make_pair(node, G[i].first[0]));
            total_vert_size += G[i].first.size();
            total_edge_size += G[i].second.size();
        }
        graph C;
        C.first.reserve(total_vert_size);
        C.second.reserve(total_edge_size);

        C.first.push_back(node);

        for(size_t i = 0; i < G.size(); i++) {
            C.first.insert(C.first.end(), G[i].first.begin(), G[i].first.end());
            C.second.insert(C.second.end(), G[i].second.begin(), G[i].second.end());
        }

        return C;
    }


    for(size_t i = 0; i < index; i++) {
        G[i].second.push_back(std::make_pair(node, G[i].first[0]));
        total_vert_size += G[i].first.size();
        total_edge_size += G[i].second.size();
    }
    for(size_t i = index + 1; i < G.size(); i++) {
        G[i].second.push_back(std::make_pair(node, G[i].first[0]));
        total_vert_size += G[i].first.size();
        total_edge_size += G[i].second.size();
    }
    total_vert_size += G[index].first.size();
    total_edge_size += G[index].second.size();
    
    graph C;
    C.first.reserve(total_vert_size);
    C.second.reserve(total_edge_size);

    for(size_t i = 0; i < G.size(); i++) {
        C.first.insert(C.first.end(), G[i].first.begin(), G[i].first.end());
        C.second.insert(C.second.end(), G[i].second.begin(), G[i].second.end());
    }

    return C;
}

graph getGraph(std::vector<int> V, std::vector<std::pair<int, int>> E, int node) {
    int index = -1;

    std::vector<graph> G = findComponents(V, E, node, &index);
    graph C = connectGraph(G, node, index);

    return C;
}

void printComponents(std::vector<graph> G) {
    std::vector<int> x;
    std::vector<std::pair<int, int>> y;
    for(size_t i = 0; i < G.size(); i++) {
        std::cout << "Graph " << i + 1 << ":\n\t";
        x = G[i].first;
        y = G[i].second;

        std::cout << "Verts: ";
        for(size_t j = 0; j < x.size(); j++) {
            std::cout << x[j] << " ";
        }
        std::cout << "\n\t";

        std::cout << "Edges: ";
        for(size_t j = 0; j < y.size(); j++) {
            std::cout << "{" << y[j].first << ", " << y[j].second << "} ";
        }
        std::cout << "\n";
    }
}

void printGraph(graph C) {
    std::cout << "Verts: ";
    for(size_t j = 0; j < C.first.size(); j++) {
        std::cout << C.first[j] << " ";
    }
    std::cout << "\n";

    std::cout << "Edges: ";
    for(size_t j = 0; j < C.second.size(); j++) {
        std::cout << "{" << C.second[j].first << ", " << C.second[j].second << "} ";
    }
    std::cout << "\n";
}

int main() {
    std::vector<int> V;
    std::vector<std::pair<int, int>> E;

    std::cout << "Enter verts: ";
    std::string line;
    std::getline(std::cin, line);

    std::istringstream sin(line);
    int num;

    while(sin >> num) {
        V.push_back(num);
    }

    std::cout << "Enter edge pairs: ";
    std::getline(std::cin, line);
    sin.clear();
    sin.str(line);
    int num1, num2;
    char b;

    while(sin >> b >> num1 >> b >> num2 >> b) {
        E.push_back(std::make_pair(num1, num2));
    }

    graph G = getGraph(V, E, 0);

    printGraph(G);
}