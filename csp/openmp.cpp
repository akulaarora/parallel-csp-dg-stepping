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

#include <omp.h>

#include "common.h"

typedef struct RelaxRequest {
    const Path_t* ali;
    const Neighbor_t* i_prime;
} RelaxRequest_t;

// This is because we do not know how many locks we need to initialize
// We will initialize them as we need them
typedef struct Lock {
    bool is_init;
    omp_lock_t lock;
} Lock_t;

// GLOBAL VARIABLES 
Bucket2D A;
Bucket3D B;

double W;
double L;
double Delta;
double Gamma;

// Parallel data structures
int num_threads = omp_get_num_threads();
// TODO buffer
std::unordered_map<int, Bucket3D> p_B; // Each processes's individual buckets
std::unordered_map<int, std::vector<RelaxRequest_t>> p_buffers; // Each processes's buffer
// Rather than creating U, we assume a thread is responsible for all nodes where
// node % total_threads == curr_thread_num. TODO we can play with this a bit.
// Locks
std::unordered_map<int, Lock_t> A_locks;
std::unordered_map<int, std::unordered_map<int, Lock_t>> B_locks;
std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, Lock_t>>> p_B_locks;
std::unordered_map<int, Lock_t> p_buffers_locks;

// HELPER FUNCTIONS
void lock(Lock_t& l) {
    if (!l.is_init) {
        omp_init_lock(&l.lock);
        l.is_init = true;
    }
    omp_set_lock(&l.lock);
}

void unlock(Lock_t& l) {
    omp_unset_lock(&l.lock);
}

// Intersection of Bjk and Uq
PathSet find_Bjk_x_Uq(int min_j, int min_k) {
    int my_thread_num = omp_get_thread_num();
    lock(B_locks[min_j][min_k]);
    PathSet ret;
    for (const Path_t& elem : B[min_j][min_k]) {
        if ((elem.path.back() % num_threads == my_thread_num)
                && (ret.count(elem) == 0)) {
            ret.insert(elem);
        }
    }
    unlock(B_locks[min_j][min_k]);
    return ret;
}


// MAIN FUNCTIONS (in paper)

void throw_req(const Path_t& ali, const Neighbor_t& i_prime) {
    // TODO need to be careful about duplicating values
    RelaxRequest_t request;
    request.ali = &ali;
    request.i_prime = &i_prime;

    int dest_thread = i_prime.node % num_threads;

    lock(p_buffers_locks[dest_thread]);
    p_buffers[dest_thread].push_back(request);
    unlock(p_buffers_locks[dest_thread]);
}

void relax(Path_t ali, Neighbor_t i_prime) {
    // TODO this is really coarse and may cause deadlocks ??
    lock(A_locks[i_prime.node]);

    Path_t new_path;
    new_path.path = ali.path;
    new_path.path.push_back(i_prime.node);
    new_path.total_cost = ali.total_cost + i_prime.cost;
    new_path.total_weight = ali.total_weight + i_prime.weight;

    bool is_dominated = false; // checks if dominated by any path in A[i']
    for (const auto& a : A[i_prime.node]) {
        if (a.total_cost < new_path.total_cost && a.total_weight <= new_path.total_weight) {
            is_dominated = true;
            break;
        }
    }
    
    if (new_path.total_weight <= W && new_path.total_cost <= L && !is_dominated) {
        A[i_prime.node].insert(new_path);
        int i = std::ceil(new_path.total_cost / Delta);
        int j = std::ceil(new_path.total_weight / Gamma);
        lock(B_locks[i][j]);
        B[i][j].insert(new_path);
        unlock(B_locks[i][j]);
    }

    auto it = A[i_prime.node].begin();
    while (it != A[i_prime.node].end()) {
        const auto& a = *it;
        if (new_path.total_cost < a.total_cost && new_path.total_weight <= a.total_weight) {
            int i = std::ceil(a.total_cost / Delta);
            int j = std::ceil(a.total_weight / Gamma);
            lock(B_locks[i][j]);
            B[i][j].erase(a);
            unlock(B_locks[i][j]);
            it = A[i_prime.node].erase(it);
        } else {
            ++it;
        }
    }

    unlock(A_locks[i_prime.node]);
}

Path_t sequential_delta_gamma_stepping(Graph& G, double inp_W, double inp_L, int start, int end,
                                       double inp_Delta, double inp_Gamma) {
    // Globalize variables
    W = inp_W;
    L = inp_L;
    Delta = inp_Delta;
    Gamma = inp_Gamma;
    
    // TODO remove this. For debugging.
    std::cout << "There are " << num_threads << " threads" << std::endl;

    Path_t initial_path;
    initial_path.path = {start};
    initial_path.total_cost = 0;
    initial_path.total_weight = 0;

    A[start].insert(initial_path);
    B[1][1].insert(initial_path);
    // Don't need to initialize locks here because our logic will initialize first time it's locked

    while (true) {
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


        PathSet R; // TODO this is temporary for heavy portion
        #pragma omp parallel default(shared) private(R)
        {
            int my_thread_num = omp_get_thread_num();
            // Line 7/8
            PathSet tmp = find_Bjk_x_Uq(min_j, min_k);
            while (!tmp.empty()) {
                // Line 9
                for (const Path_t& elem : tmp) {
                    if (R.count(elem) == 0) {
                        R.insert(elem);
                    }
                }
                // TODO this could be optimized (put a barrier perhaps and clear B all at once)
                // Line 10
                lock(B_locks[min_j][min_k]);
                for (const Path_t& elem : tmp) {
                    B[min_j][min_k].erase(elem);
                }
                unlock(B_locks[min_j][min_k]);

                // Line 11-12
                for (const Path_t& ali : tmp) {
                    int i = ali.path.back();
                    for (const Neighbor_t& i_prime : G.neighbor(i)) {
                        if (i_prime.cost < Delta && i_prime.weight < Gamma) {
                            throw_req(ali, i_prime);
                        }
                    }
                }

                // Line 13-14
                lock(p_buffers_locks[my_thread_num]);
                for (const RelaxRequest_t& req : p_buffers[my_thread_num]) {
                    relax(*req.ali, *req.i_prime);
                }
                unlock(p_buffers_locks[my_thread_num]);
            }

        }

        // TODO WE HAVE NOT PARALLELIZED THIS YET
        for (const auto& ali : R) {
            int i = ali.path.back();
            for (const Neighbor_t& i_prime : G.neighbor(i)) {
                if (i_prime.cost >= Delta || i_prime.weight >= Gamma) {
                    relax(ali, i_prime);
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

