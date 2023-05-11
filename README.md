# Parallelizing Constrained Shortest Path Using (∆, Γ)-Stepping Algorithm

In this paper, we implement a serial implementation and a shared-memory parallel implementation of
the (∆, Γ)-stepping algorithm to solve the NP-hard Constrained Shortest Path (CSP) problem. The
CSP problem says that for a directed graph with non-negative costs and weights for all edges in the
graph, find the shortest path (by cost) such that the total weight of the path is less than a given positive
integer W. CSP is an NP-hard problem and a harder variant of Single-Source Shortest Path (solved
using the ∆-Stepping parallel algorithm). The (∆, Γ)-stepping algorithm is a theoretical algorithm
proposed by Bahreini et al. (2022) that has never been implemented before. Our serial and parallel
implementations are both novel, having never been attempted before. Our serial implementation
is quite effective by solving graphs up to 10,000 edges and nodes, which a naive algorithm cannot
complete in feasible time. Our parallel implementation demonstrates a major improvement over the
serial and naive algorithms by solving graphs up to 1 million edges and nodes, which neither the serial
nor naive algorithm can achieve in a feasible time. On a graph with 10,000 edges and nodes, our
parallel algorithm achieves a 70x speedup over the serial implementation. We additionally discover an
issue and an optimization in the theoretical algorithm, which we discuss and fix in our implementation.
We also discuss some of the nuances of this algorithm, in particular the implementation details, the
performance and the performance’s analysis, and the speed profile of the parallel algorithm. Lastly,
we provide some intuition for why we believe our parallel implementation will perform better than a
distributed-memory implementation and provide some areas for further research, which we hope can
build off our work to continue to push the boundary of parallel graph algorithms.

## Code

dg.py provides a Python implementation of the serial algorithm.

The csp folder provides a C++ implementation of the algorithm as well as a pparallel implementation.
