#include <unordered_map>
#include <vector>

#include "common.h"

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

const std::vector<Neighbor_t>& Graph::neighbor(int node) {
    static const std::vector<Neighbor_t> empty_list;
    if (neighbors.find(node) != neighbors.end()) {
        return neighbors[node];
    } else {
        return empty_list;
    }
}