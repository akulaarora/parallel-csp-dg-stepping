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
    Path_t ali;
    Neighbor_t i_prime;
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
int num_threads = 64;
// TODO buffer
std::unordered_map<int, std::vector<RelaxRequest_t>> p_buffers; // Each processes's buffer
// Rather than creating U, we assume a thread is responsible for all nodes where
// node % total_threads == curr_thread_num. TODO we can play with this a bit.
// Locks
// std::unordered_map<int, std::unordered_map<int, Lock_t>> B_locks; // These are dynamically initialized
Lock_t B_lock;
std::unordered_map<int, Lock_t> p_buffers_locks; // Statically initialized
Lock_t init_lock; // Lock for initializing. We need this to ensure two threads do not initialize the same lock.

// HELPER FUNCTIONS
void lock(Lock_t& l) {
    // if (!l.is_init) {
    //     omp_set_lock(&init_lock.lock);
    //     // Double checking to ensure it wasn't initialized by another thread
    //     // while we were checking/waiting for the lock
    //     if (!l.is_init) {
    //         omp_init_lock(&l.lock);
    //         l.is_init = true;
    //     }
    //     omp_unset_lock(&init_lock.lock);
    // }
    omp_set_lock(&l.lock);
}

void unlock(Lock_t& l) {
    omp_unset_lock(&l.lock);
}

// Intersection of Bjk and Uq
PathSet find_Bjk_x_Uq(int min_j, int min_k) {
    int my_thread_num = omp_get_thread_num();
    PathSet ret;
    for (const Path_t& elem : B[min_j][min_k]) {
        if ((elem.path.back() % num_threads == my_thread_num)
                && (ret.count(elem) == 0)) {
            ret.insert(elem);
        }
    }
    return ret;
}

// MAIN FUNCTIONS (in paper)

void throw_req(const Path_t& ali, const Neighbor_t& i_prime) {
    // TODO need to be careful about duplicating values
    RelaxRequest_t request;
    request.ali = ali;
    request.i_prime = i_prime;

    int dest_thread = i_prime.node % num_threads;

    lock(p_buffers_locks[dest_thread]);
    p_buffers[dest_thread].push_back(request);
    unlock(p_buffers_locks[dest_thread]);
}

// TODO can make this pass by reference ?
void relax(Path_t ali, Neighbor_t i_prime) {
    // TODO can make B locking finer to per bucket

    Path_t new_path;
    new_path.path = ali.path;
    new_path.path.push_back(i_prime.node);
    new_path.total_cost = ali.total_cost + i_prime.cost;
    new_path.total_weight = ali.total_weight + i_prime.weight;

    // Print out new path
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
        lock(B_lock);
        B[i][j].insert(new_path);
        unlock(B_lock);
    }

    auto it = A[i_prime.node].begin();
    while (it != A[i_prime.node].end()) {
        const auto& a = *it;
        if (new_path.total_cost < a.total_cost && new_path.total_weight <= a.total_weight) {
            int i = std::ceil(a.total_cost / Delta);
            int j = std::ceil(a.total_weight / Gamma);
            lock(B_lock);
            B[i][j].erase(a);
            unlock(B_lock);
            it = A[i_prime.node].erase(it);
        } else {
            ++it;
        }
    }
}

Path_t sequential_delta_gamma_stepping(Graph& G, double inp_W, double inp_L, int start, int end,
                                       double inp_Delta, double inp_Gamma) {
    // Globalize variables
    W = inp_W;
    L = inp_L;
    Delta = inp_Delta;
    Gamma = inp_Gamma;
    
    Path_t initial_path;
    initial_path.path = {start};
    initial_path.total_cost = 0;
    initial_path.total_weight = 0;

    for (int i = 0; i < G.total_nodes; ++i) {
        // We don't need locks for this bc each thread will only access its own A[i']
        A[i] = PathSet();
    }

    A[start].insert(initial_path);
    B[1][1].insert(initial_path);

    // Initialize static locks
    // omp_init_lock(&init_lock.lock);
    // init_lock.is_init = true;

    for (int i = 0; i < num_threads; ++i) {
        omp_init_lock(&p_buffers_locks[i].lock);
        p_buffers_locks[i].is_init = true;
        p_buffers[i] = std::vector<RelaxRequest_t>(); // Add this line to initialize p_buffers for each thread
    }

    omp_init_lock(&B_lock.lock);
    B_lock.is_init = true;

    // Line 3
    while (true) {
        // Line 4
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

        /*
         * IMPORTANT: There is an issue with the paper that lines 13-14 are 
         * in the while loop conditioned on B[j][k] x Uq =/= 0 for processor q.
         * Processor q may not have any elements in B[j][k] x Uq, but it may
         * still have elements added to its buffer from other processors.
         * TODO we should write out the correct paper algorithm.
         * 
         * My solution is to change line 7 to while B[j][k] is not empty
         * (all processors run in this while loop), then proceed thru
         * lines 8-12 if B[j][k] x Uq is not empty for processor q.
         * All processors will then run lines 13-14.
         * 
         * Lines marked as "NEW" are the corrected algo lines.
         */

        // Line 5
        #pragma omp parallel default(shared) num_threads(64)
        {
            // num_threads = omp_get_num_threads();
            int my_thread_num = omp_get_thread_num();

            // Line 6
            PathSet R;

            // NEW Line 7
            while (!B[min_j][min_k].empty()) {
                // Line 8
                PathSet tmp;
                tmp = find_Bjk_x_Uq(min_j, min_k);

                // Line 9
                for (const Path_t& elem : tmp) {
                    if (R.count(elem) == 0) {
                        R.insert(elem);
                    }
                }

                #pragma omp barrier

                // TODO this could be optimized (put a barrier perhaps and clear B all at once)
                // Line 10
                // TODO see if we can get rid of this and add more synchronization (what paper describes)
                // NOTE: Removed locks in find_Bjk_x_Uq bc there's no modification of B while that runs
                // with this implementation
                #pragma omp master
                {
                    // No need for lock bc only one thread is here
                    B[min_j][min_k].clear();
                }
            
                // Line 11-12
                for (const Path_t& ali : tmp) {
                    int i = ali.path.back();
                    for (const Neighbor_t& i_prime : G.neighbor(i)) {
                        if (i_prime.cost < Delta && i_prime.weight < Gamma) {
                            throw_req(ali, i_prime);
                        }
                    }
                }

                // TODO there may be a way that's closer to the essence of the paper
                // that allows for synchronous throwing and relaxing
                #pragma omp barrier

                // Line 13-14
                lock(p_buffers_locks[my_thread_num]);
                for (const RelaxRequest_t& req : p_buffers.at(my_thread_num)) {
                    relax(req.ali, req.i_prime);
                }
                // Clear the buffer
                p_buffers[my_thread_num].clear();
                unlock(p_buffers_locks[my_thread_num]);

                #pragma omp barrier
            }

            for (const Path_t& ali : R) {
                int i = ali.path.back();
                for (const Neighbor_t& i_prime : G.neighbor(i)) {
                    if (i_prime.cost >= Delta || i_prime.weight >= Gamma) {
                        throw_req(ali, i_prime);
                    }
                }
            }

            #pragma omp barrier

            // Line 13-14
            lock(p_buffers_locks[my_thread_num]);
            for (const RelaxRequest_t& req : p_buffers.at(my_thread_num)) {
                relax(req.ali, req.i_prime);
            }
            // Clear the buffer
            p_buffers[my_thread_num].clear();
            unlock(p_buffers_locks[my_thread_num]);
                
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

