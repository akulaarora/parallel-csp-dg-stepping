#include "common.h"
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <limits>

void relax(Bucket3D& B, 
           Bucket2D& A, 
           Path_t ali, Neighbor_t i_prime, double W, double L, double Delta, double Gamma) {
    Path_t new_path;
    new_path.path = ali.path;
    new_path.path.push_back(i_prime.node);
    new_path.total_cost = ali.total_cost + i_prime.cost;
    new_path.total_weight = ali.total_weight + i_prime.weight;

    bool is_dominated = false;
    for (const auto& a : A[i_prime.node]) {
        if (a.total_cost < new_path.total_cost && a.total_weight <= new_path.total_weight) {
            is_dominated = true;
            break;
        }
    }

    if (new_path.total_weight <= W && new_path.total_cost <= L && !is_dominated) {
        A[i_prime.node].insert(new_path);
        B[std::ceil(new_path.total_cost / L)][std::ceil(new_path.total_weight / W)].insert(new_path);
    }

    for (const auto& a : A[i_prime.node]) {
        if (new_path.total_cost < a.total_cost && new_path.total_weight <= a.total_weight) {
            A[i_prime.node].erase(a);
            B[std::ceil(a.total_cost / L)][std::ceil(a.total_weight / W)].erase(a);
        }
    }
}

Path_t sequential_delta_gamma_stepping(Graph& G, double W, double L, int start, int end, double Delta, double Gamma) {
    Bucket2D A;
    Bucket3D B;

    Path_t initial_path;
    initial_path.path = {start};
    initial_path.total_cost = 0;
    initial_path.total_weight = 0;

    A[start].insert(initial_path);
    B[1][1].insert(initial_path);

    while (true) {
        int min_j = -1;
        int min_k = -1;

        for (int j = 1; j <= std::ceil(L / Delta); ++j) {
            for (int k = 1; k <= std::ceil(W / Gamma); ++k) {
                if (!B[j][k].empty() && (min_j == -1 || min_k == -1 || (j < min_j) || (j == min_j && k < min_k))) {
                    min_j = j;
                    min_k = k;
                }
            }
        }

        if (min_j == -1 && min_k == -1) {
            break;
        }

        PathSet R;
        while (!B[min_j][min_k].empty()) {
            for (const auto& elem : B[min_j][min_k]) {
                R.insert(elem);
            }
            B[min_j][min_k].clear();

            for (const auto& ali : R) {
                int i = ali.path.back();
                for (const auto& i_prime : G.neighbor(i)) {
                    if (i_prime.cost < Delta && i_prime.weight < Gamma) {
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
        return Path_t{{}, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()};
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

int main() {
    return 1;
}