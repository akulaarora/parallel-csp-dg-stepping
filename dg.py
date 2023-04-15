import math
import networkx as nx

class Algorithm:
    def __init__(self):
        self.distances = {}
        self.delta = 5
        self.gamma = 5
        self.property_map = {}
        self.source_vertex = 2
        self.infinity = float("inf")
        self.B = {}
        self.W = 0
        self.L = 0

        self.A = {}

    def relax(self, a, i_dest, g):
        '''
        a^l_i is a label on node i
        i_dest is the destination node
        '''
        # TODO
        c_new = a[1] + g.get_edge_data(a[0], i_dest)['cost']
        w_new = a[2] + g.get_edge_data(a[0], i_dest)['weight']
        a_new = (i_dest, )

        if w_new <= self.W and c_new <= self.L:
            updated = False
            to_remove = []
            for a_prime in self.B[w]:
                if x > a_prime:
                    updated = True
                elif x < a_prime:
                    to_remove.append(a_prime)

            if updated:
                self.B[w].add(x)

                for a in to_remove:
                    self.B[w].remove(a)

    def find_requests(self, vertices, kind, g):
        """
        returns a dictionary of neighboring edges with their weights but according to the kind i.e. light or heavy

        :param vertices:
        :param kind:
        :param g:
        :return:
        """

        tmp = {}
        # print("vertices=", vertices, "kind=", kind)
        for u in vertices:
            for v in g.neighbors(u):
                # print(u, self.property_map[u], g.get_edge_data(u, v)['weight'])
                edge_weight = self.property_map[u] + g.get_edge_data(u, v)['weight']
                if kind == 'light':
                    if g.get_edge_data(u, v)['weight'] <= self.delta:
                        if v in tmp:
                            if edge_weight < tmp[v]:
                                tmp[v] = edge_weight
                        else:
                            tmp[v] = edge_weight
                elif kind == 'heavy':
                    if g.get_edge_data(u, v)['weight'] > self.delta:
                        if v in tmp:
                            if edge_weight < tmp[v]:
                                tmp[v] = edge_weight
                        else:
                            tmp[v] = edge_weight
                else:
                    return "Error: No such kind of edges " + kind
        # print("tmp=", tmp)
        return tmp

    def find_requests(self, u, kind, g):
        # u is the vertex
        tmp = {}
        for v in g.neighbors(u):
            edge_cost = g.get_edge_data(u, v)['cost']
            edge_weight = g.get_edge_data(u, v)['weight']
            edge_value = (v, self.property_map[u][0] + edge_cost, self.property_map[u][1] + edge_weight)

            # TODO not sure about light being below edge cost/weight
            if kind == 'light' and edge_cost <= self.delta and edge_weight <= self.gamma:
                if v in tmp:
                    if edge_value < tmp[v]:
                        tmp[v] = edge_value
                else:
                    tmp[v] = edge_value
            elif kind == 'heavy':
                if v in tmp:
                    if edge_value < tmp[v]:
                        tmp[v] = edge_value
                else:
                    tmp[v] = edge_value
            else:
                return "Error: No such kind of edges " + kind
        return tmp

    def relax_requests(self, request):
        for key, value in request.items():
            self.relax(key, value)

    def min_lex(self, B):
        for i in range(math.ceil(self.W / self.delta)):
            for j in range(math.ceil(self.L / self.gamma)):
                if B[i][j]:
                    return (i, j)
        return None

    def delta_gamma_stepping(self, g, W, L):
        self.W = W
        self.L = L

        for node in g.nodes():
            self.A[node] = set()
        self.A[self.source_vertex].add((self.source_vertex, 0, 0)) # TODO

        for i in range(math.ceil(W / self.delta)):
            self.B[i] = {}
            for j in range(math.ceil(L / self.gamma)):
                self.B[i][j] = set()
        self.B[1][1].add((self.source_vertex, 0, 0))

        while True:
            min_lex = self.min_lex(self.B)
            if min_lex is None:
                break

            R = {}
            
            Bjk = self.B[min_lex[0]][min_lex[1]]
            while len(self.B[min_lex[0]][min_lex[1]]) > 0:
                R = {**R, **self.B[min_lex[0][min_lex[1]]]}
                tmp = dict(Bjk)
                self.B[min_lex[0]][min_lex[1]] = {}
                
                for a in tmp:
                    self.relax(a[0], (a[1], a[2]))

                for 

                u, c, w = self.B[min_lex[0]][min_lex[1]].pop()
                R.add(u)
                for v in g.neighbors(u):
                    edge_cost = g.get_edge_data(u, v)['cost']
                    edge_weight = g.get_edge_data(u, v)['weight']
                    if w + edge_weight <= W and c + edge_cost <= L:
                        self.B[math.ceil((w + edge_weight) / self.delta)][math.ceil((c + edge_cost) / self.gamma)].add((v, c + edge_cost, w + edge_weight))

        while any(self.B[node] for node in g.nodes()):
            r = set()
            while True:
                b = {node: min(self.B[node]) for node in g.nodes() if self.B[node]}
                if not b:
                    break

                i = min(b, key=b.get)
                x = b[i]

                self.B[i].remove(x)
                r.add(i)
                req = self.find_requests({i}, 'light', g)
                self.relax_requests(req)

            req = self.find_requests(r, 'heavy', g)
            self.relax_requests(req)


    def validate(self, g, target_algorithm):
        # TODO
        pass
        # self.property_map = {k: v for k, v in self.property_map.items() if v != (self.infinity, self.infinity)}
        # valid = True

        # for node in g.nodes():
        #     if target_algorithm(self.property_map[node]) != target_algorithm(self.distances[node]):
        #         print("Error: The algorithm is faulty!")
        #         print("vertex ", node, " value in ground truth is ", self.distances[node], " and value in delta-gamma stepping is ", self.property

