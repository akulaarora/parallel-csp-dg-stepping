import math
from collections import defaultdict
import numpy as np
import networkx as nx
from networkx.algorithms.simple_paths import all_simple_paths

def get_route_edge_attributes(
    G, route, attribute=None, minimize_key="cost", retrieve_default=None
):
    attribute_values = []
    for u, v in zip(route[:-1], route[1:]):
        # if there are parallel edges between two nodes, select the one with the
        # lowest value of minimize_key
        data = min(G.get_edge_data(u, v).values(), key=lambda x: x[minimize_key])
        if attribute is None:
            attribute_value = data
        elif retrieve_default is not None:
            attribute_value = data.get(attribute, retrieve_default(u, v))
        else:
            attribute_value = data[attribute]
        attribute_values.append(attribute_value)
    return attribute_values

class Graph:
    # ["A", "B", "C", "D", "E", "F", "G"]
    # [1, 2, 3, 4, 5, 6, 7] # vertices are numbers from 1 to n
    # neighbors (node: 1 to n) -> [(node, cost, weight), ...] 
    '''
    array or a vector ->
    struct Neighbor {
        int node;
        double cost;
        double weight;
    };
    
    Graph G:
        public:
            vector<int> neighbors;
            neighbor(int) -> vector<Neighbor>;
    '''
    def __init__(self, MG):
        self.MG = MG
        self.nodes = {i+1: node for i, node in enumerate(MG.nodes)}
        print(self.nodes)
        self.inv_nodes = {v: k for k, v in self.nodes.items()}
 
    def neighbors(self, node):
        node_inv = self.nodes[node]
        neighbors_inv = list(self.MG.neighbors(self.nodes[node]))
        neighbors = [self.inv_nodes[neighbor] for neighbor in neighbors_inv]
        ret = []
        for i, neighbor_inv in enumerate(neighbors_inv):
            vals = get_route_edge_attributes(self.MG, [node_inv, neighbor_inv])[0]
            ret.append((neighbors[i], vals['cost'], vals['weight']))
        return ret

class Algorithm:
    def __init__(self, Delta, Gamma):
        self.Delta = Delta
        self.Gamma = Gamma

    def relax(self, B, A, ali, i_prime, W, L):
        '''
        Relaxes an edge that we can now include in the CSP calculation
        and updates paths accordingly
        '''
        pi, c, w = ali
        pi_new, c_new, w_new = pi + (i_prime[0],), c + i_prime[1], w + i_prime[2] # building on the path

        is_dominated=False
        for a in A[i_prime[0]]:
            if a[1] < c_new and a[2] <= w_new:
                is_dominated = True
                break
                
        if w_new <= W and c_new <= L and not is_dominated:
            A[i_prime[0]].add((pi_new, c_new, w_new))
            B[math.ceil(c_new / L)][math.ceil(w_new / W)].add((pi_new, c_new, w_new))

        
        for a in A[i_prime[0]]:
            if c_new < a[1] and w_new <= a[2]:
                A[i_prime[0]].remove(a)
                B[math.ceil[a[1] / L][math.ceil[a[2] / W]]].remove(a)

    def sequential_delta_gamma_stepping(self, MG, W, L, start, end):
        # W is the constraint on weight and L is the constraint on cost
        # we shouldn't use L (set it to max double) since it's not part of CSP
        G = Graph(MG)
        start = G.inv_nodes[start]
        end = G.inv_nodes[end]
        A = defaultdict(set) # initialize from 1 to n
        B = defaultdict(lambda: defaultdict(set)) # vector<vector<Hashset>> 
        # array instead of vector
        # main vector initialize from 0 to (L+1)/self.Delta, nested vector 0 to (W+1)/self.Gamma

        '''
        # this or a tuple of (LinkedList/vector<int>, float, float)
        struct Path {
            LinkedList<int> path; // list of all nodes in the path
            float total_cost;
            float total_weight;
        } PATH_T;
        '''
        A[1].add(((start,), 0, 0))
        B[1][1].add(((start,), 0, 0))

        # Check if any bucket is non-empty (B[j][k] for all j, k)
        while any(B[j][k] for j in range(1, math.ceil(L/self.Delta) + 1) for k in range(1, math.ceil(W/self.Gamma) + 1)):
            # Getting the minimum non-empty bucket (lexicographic min by j then k)
            j, k = min((j, k) for j in range(1, math.ceil(L/self.Delta) + 1) for k in range(1, math.ceil(W/self.Gamma) + 1) if B[j][k])

            R = set() # hashset
            while B[j][k]:
                R |= B[j][k]
                tmp = B[j][k]
                B[j][k] = set()

                for ali in tmp:
                    i = ali[0][-1] # the last node in the path
                    for i_prime in G.neighbors(i):
                        if i_prime[1] < self.Delta and i_prime[2] < self.Gamma:
                            self.relax(B, A, ali, i_prime, W, L)

            for ali in R:
                i = ali[0][-1] # the last node in the path
                for i_prime in G.neighbors(i):
                    if i_prime[1] >= self.Delta or i_prime[2] >= self.Gamma:
                        self.relax(B, A, ali, i_prime, W, L)

        if not A[end]:
            # return struct path or tuple of (LinkedList<int>, double, double)
            return [], float('inf'), float('inf')
        else:
            # something that is not valid or throw an exception
            # return struct path or tuple of (null, -1, -1);
            # return null
            return min(A[end], key=lambda x: x[1])


def validate(MG):
    # something like this in cpp or just do in python and compare the results
    shortest_paths_constrained = []
    for path in all_simple_paths(MG, "A", "G"):
        if sum(get_route_edge_attributes(MG, path, 'weight')) < 300:
            path_cost = sum(get_route_edge_attributes(MG, path, 'cost'))
            result = (path, path_cost)
            if not shortest_paths_constrained:
                shortest_paths_constrained.append(result)
            elif path_cost == shortest_paths_constrained[0][1]:
                shortest_paths_constrained.append(result)
            elif path_cost < shortest_paths_constrained[0][1]:
                shortest_paths_constrained = [result]

    return shortest_paths_constrained

def main():
    # 1. Define test network 
    MG = nx.MultiDiGraph()
    MG.add_edges_from([("B", "A", {"cost": 0.8}), ("A", "B", {"cost": 1.}), ("D", "G", {"cost": 3.5}),
                    ("B", "D", {"cost": 20.8}), ("A", "C", {"cost": 9.7}), ("D", "C", {"cost": 0.3}),
                    ("B", "E", {"cost": 4.8}), ("D", "E", {"cost": 0.05}), ("C", "E", {"cost": 0.1}),          
                    ("E", "C", {"cost": 0.7}), ("E", "F", {"cost": 0.4}), ("E", "G", {"cost": 15.}),           
                    ("F", "C", {"cost": 0.9}), ("F", "D", {"cost": 4.}),                       
                    ])
    attrs = {'B': {"x": -20., "y": 60.}, 'A': {"x": 28., "y":55.},
            'C': {"x": -12., "y": 40.}, 'D': {"x": 40., "y":45.},
            'E': {"x": 8., "y": 35.}, 'F': {"x": -8., "y":15.},    
            'G': {"x": 21., "y":5.},    
            }

    for num, (k,v) in enumerate(attrs.items()):
        attrs[k]={**v,
                }  
    nx.set_node_attributes(MG, attrs)

    rng = np.random.default_rng(seed=42)
    random_weight = list(rng.uniform(low=0, high=100, size=len(MG.edges)).round(0))
    attrs={}
    for num, edge in enumerate(MG.edges):
        attrs[edge]={'weight': random_weight[num]}
    nx.set_edge_attributes(MG, attrs)

    print(validate(MG))

    x = Algorithm(10, 10)
    print(x.sequential_delta_gamma_stepping(MG, 300, 99999, "A", "G"))

main()