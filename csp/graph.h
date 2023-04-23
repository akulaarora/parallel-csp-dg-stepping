#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <unordered_map>
#include <vector>

#include "common.h"

class Graph {
    private:
        std::unordered_map<int, std::vector<Neighbor_t>> neighbors;
    public:
        Graph(int num_nodes, const std::vector<Edge_t>& edges);
        int total_nodes;
        std::vector<Neighbor_t> neighbor(int node);
};

Graph::Graph(int num_nodes, const std::vector<Edge_t>& edges) {
    total_nodes = num_nodes;
    for (int i = 0; i != edges.size(); i++) {
        if (neighbors.find(edges[i].src) == neighbors.end()) {
            neighbors[edges[i].src] = std::vector<Neighbor_t>();
        }
        Neighbor_t neighbor = {};
        neighbor.node = edges[i].dest;
        neighbor.cost = edges[i].cost;
        neighbor.weight = edges[i].weight;
        neighbors[edges[i].src].push_back(neighbor);
    }
}

std::vector<Neighbor_t> Graph::neighbor(int node) {
    if (neighbors.find(node) != neighbors.end()) {
        return neighbors[node];
    } else {
        return std::vector<Neighbor_t>();
    }
}

#endif
