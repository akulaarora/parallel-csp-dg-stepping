#ifndef __COMMON_H__
#define __COMMON_H__

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstddef>

typedef struct Neighbor {
    int node;
    double cost;
    double weight;
} Neighbor_t;

typedef struct Edge {
    int src;
    int dest;
    double cost;
    double weight;
} Edge_t;

typedef struct Path {
    std::vector<int> path; // list of all nodes in the path
    double total_cost;
    double total_weight;
} Path_t;

struct PathHash {
    size_t operator()(const Path_t& p) const {
        std::vector<int> v = p.path;
        std::hash<int> hasher;
        size_t seed = 0;
        for (int i : v) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }
        return seed;
    }
};

struct PathEqual {
    bool operator()(const Path_t& lhs, const Path_t& rhs) const {
        return lhs.path == rhs.path && lhs.total_cost == rhs.total_cost && lhs.total_weight == rhs.total_weight;
    }
};

using PathSet = std::unordered_set<Path_t, PathHash, PathEqual>;
using Bucket2D = std::unordered_map<int, PathSet>;
using Bucket3D = std::unordered_map<int, Bucket2D>;

#endif
