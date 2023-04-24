#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <time.h>
#include <unordered_map>
#include <vector>

#include "common.h"

void relax(Bucket3D& B, Bucket2D& A, Path_t ali, Neighbor_t i_prime, double W, double L,
           double Delta, double Gamma) {
#ifdef DEBUG
    std::cout << "RELAXING STARTS" << std::endl;
#endif

    Path_t new_path;
    new_path.path = ali.path;
    new_path.path.push_back(i_prime.node);
    new_path.total_cost = ali.total_cost + i_prime.cost;
    new_path.total_weight = ali.total_weight + i_prime.weight;

#ifdef DEBUG
    std::cout << "New Path: ";
    for (const auto& node : new_path.path) {
        std::cout << node << " ";
    }
    std::cout << std::endl;
    std::cout << "New Path Cost: " << new_path.total_cost << std::endl;
    std::cout << "New Path Weight: " << new_path.total_weight << std::endl;
#endif

    bool is_dominated = false; // checks if dominated by any path in A[i']
    for (const auto& a : A[i_prime.node]) {
        if (a.total_cost < new_path.total_cost && a.total_weight <= new_path.total_weight) {
            is_dominated = true;
            break;
        }
    }
    
#ifdef DEBUG
    std::cout << "is_dominated: " << is_dominated << std::endl;
#endif

    if (new_path.total_weight <= W && new_path.total_cost <= L && !is_dominated) {
#ifdef DEBUG
        std::cout << "Inserting into A and B" << std::endl;
#endif

        A[i_prime.node].insert(new_path);
#ifdef DEBUG
        std::cout << std::ceil(new_path.total_cost / Delta) << std::endl;
        std::cout << std::ceil(new_path.total_weight / Gamma) << std::endl;
#endif
        B[std::ceil(new_path.total_cost / Delta)][std::ceil(new_path.total_weight / Gamma)].insert(
            new_path);
    }

    auto it = A[i_prime.node].begin();
    while (it != A[i_prime.node].end()) {
        const auto& a = *it;
        if (new_path.total_cost < a.total_cost && new_path.total_weight <= a.total_weight) {
#ifdef DEBUG
            std::cout << "Erasing Path from A[i_prime.node]: ";
            for (const auto& node : a.path) {
                std::cout << node << " ";
            }
            std::cout << std::endl;
            std::cout << "Erased Path Cost: " << a.total_cost << std::endl;
            std::cout << "Erased Path Weight: " << a.total_weight << std::endl;
#endif
            B[std::ceil(a.total_cost / Delta)][std::ceil(a.total_weight / Gamma)].erase(a);
            it = A[i_prime.node].erase(it);
        } else {
            ++it;
        }
    }
}

Path_t sequential_delta_gamma_stepping(Graph& G, double W, double L, int start, int end,
                                       double Delta, double Gamma) {
    Bucket2D A;
    Bucket3D B;

    Path_t initial_path;
    initial_path.path = {start};
    initial_path.total_cost = 0;
    initial_path.total_weight = 0;

    A[start].insert(initial_path);
    B[1][1].insert(initial_path);

    while (true) {
        // std::cout << "Run" << std::endl;
        int min_j = -1;
        int min_k = -1;

        for (auto it = B.begin(); it != B.end(); ++it) {
            int j = it->first;
            for (auto it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
                int k = it2->first;
                if (!it2->second.empty() &&
                    (min_j == -1 || min_k == -1 || (j < min_j) || (j == min_j && k < min_k))) {
                    min_j = j;
                    min_k = k;
                }
            }
        }

        if (min_j == -1 && min_k == -1) {
            break;
        }

#ifdef DEBUG
        std::cout << "Going thru bucket " << min_j << ", " << min_k << std::endl;
#endif

        PathSet R;
        while (!B[min_j][min_k].empty()) {
            for (const auto& elem : B[min_j][min_k]) {
                if (R.count(elem) == 0) {
                    R.insert(elem);
                }
            }
            PathSet tmp(B[min_j][min_k]);
            B[min_j][min_k].clear();

#ifdef DEBUG
            std::cout << "R size: " << R.size() << std::endl;
            std::cout << "tmp size: " << tmp.size() << std::endl;
#endif
            for (const auto& ali : tmp) {
#ifdef DEBUG
                std::cout << "Path iterated over: ";
                for (const auto& node : ali.path) {
                    std::cout << node << " ";
                }
                std::cout << std::endl;
#endif

                int i = ali.path.back();
#ifdef DEBUG
                std::cout << "i--going thru neighbours: " << i << std::endl;
#endif
                for (const auto& i_prime : G.neighbor(i)) {
                    if (i_prime.cost < Delta && i_prime.weight < Gamma) {
#ifdef DEBUG
                        std::cout << "Relaxing " << i_prime.node << std::endl;
#endif
                        relax(B, A, ali, i_prime, W, L, Delta, Gamma);
                    }
                }
            }
        }

        for (const auto& ali : R) {
            int i = ali.path.back();
            for (const auto& i_prime : G.neighbor(i)) {
                if (i_prime.cost >= Delta || i_prime.weight >= Gamma) {
                    relax(B, A, ali, i_prime, W, L, Delta, Gamma);
                }
            }
        }
    }

    if (A[end].empty()) {
        return Path_t{
            {}, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
    } else {
        Path_t min_path = *A[end].begin();
        for (const auto& path : A[end]) {
            if (path.total_cost < min_path.total_cost) {
                min_path = path;
            }
        }
        return min_path;
    }
}
