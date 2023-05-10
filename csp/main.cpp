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

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Missing input file!\n";
        return EXIT_FAILURE;
    }

    bool weightsGiven = false;
    if (argc == 3) {
        std::string weights(argv[2]);
        if (weights == "--nogen") {
            weightsGiven = true;
        }
    }

    std::string input_file(argv[1]);
    std::ifstream input;
    std::vector<std::string> lines;
    input.open(input_file);
    std::string line;
    while (std::getline(input, line)) {
        lines.push_back(line);
    }
    input.close();

    double low = 0;
    double high = 100;
    const long max_rand = 1000000L;
    srandom(time(NULL));
    std::vector<Edge_t> edges;
    int num_nodes;
    int num_edges;

    std::ofstream output("full.txt");

    for (int i = 0; i != lines.size(); ++i) {
        std::istringstream stream(lines[i]);
        std::vector<std::string> temp;
        std::string curr;

        while (stream >> curr) {
            temp.push_back(curr);
        }

        if (i == 0) {
            num_nodes = std::stoi(temp[0]);
            num_edges = std::stoi(temp[1]);
            output << lines[i] << "\n";
        } else {
            double cost = 0.0;
            double weight = 0.0;
            if (!weightsGiven) {
                cost = low + (high - low) * (random() % max_rand) / max_rand;
                // std::cout << cost << std::endl;
                weight = low + (high - low) * (random() % max_rand) / max_rand;
                output << lines[i] << "\t" << cost << "\t" << weight << "\n";
            } else {
                cost = std::stod(temp[2]);
                weight = std::stod(temp[3]);
                output << lines[i] << "\n";
            }

            Edge_t edge = {};
            edge.src = std::stoi(temp[0]);
            edge.dest = std::stoi(temp[1]);
            edge.cost = cost;
            edge.weight = weight;
            edges.push_back(edge);
        }
    }

    output.close();

    std::cout << edges.size() << "\n";
    
    Graph graph(num_nodes, edges);

    // std::vector<Neighbor_t> node_0_neighbors = graph.neighbor(0);
    // for (const auto& neighbor : node_0_neighbors) {
    //     std::cout << "Node 0 is connected to node " << neighbor.node << " with cost " << neighbor.cost << " and weight " << neighbor.weight << "\n";
    // }

    std::cout << "Sequential Delta Gamma Stepping\n";

    auto start = std::chrono::high_resolution_clock::now();
    
    // We don't want to add L constraint since it's not part of CSP
    Path_t result = sequential_delta_gamma_stepping(graph, 300, std::numeric_limits<double>::max(), 1, 2, 10, 10);

    auto end = std::chrono::high_resolution_clock::now();

    // COMMENT OUT WHEN TIMING
    std::cout << "Cost: " << result.total_cost << std::endl;
    std::cout << "Weight: " << result.total_weight << std::endl;
    // Print path
    for (const auto& node : result.path) {
        std::cout << node << " ";
    }
    std::cout << " " << std::endl;

    std::chrono::duration<double> total = end - start;
    std::cout << "Runtime: " << total.count() << std::endl;

    return 0;
}