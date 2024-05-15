'''
Graph represented as adjacency list using a dictionary
keys are vertices
values of the dictionary represent the list of adjacent vertices of the key node
'''

class MyGraph:
    def __init__(self, g = None):
        ''' takes dictionary to fill the graph as input; default is empty dictionary '''
        if g is None:
            g = {}
        self.graph = g

    def print_graph(self):
        ''' Prints the content of the graph as adjacency list '''
        for k, v in self.graph.items():
            print (k, " -> ", v)

    ## get basic info
    def get_nodes(self):
        ''' Returns list of nodes in the graph '''
        return list(self.graph.keys())


    def get_edges(self):
        ''' Returns edges in the graph as a list of tuples (origin, destination) '''
        return [(k, v) for k in self.graph.keys() for v in self.graph[k]]

    def size(self):
        ''' Returns size of the graph : number of nodes, number of edges '''
        return len(self.graph), len(self.get_edges())

    ## add nodes and edges
    def add_vertex(self, v):
        ''' Add a vertex to the graph; tests if vertex exists not adding if it does '''
        if v not in self.graph:
            self.graph[v] = []


    def add_edge(self, o, d):
        ''' Add edge to the graph; if vertices do not exist, they are added to the graph '''
        self.add_vertex(o)
        self.add_vertex(d)
        if d not in self.graph[o]:
            self.graph[o].append(d)

    ## successors, predecessors, adjacent nodes

    def get_successors(self, v):
        # list() needed to avoid list being overwritten of result of the function is used
        return list(self.graph[v])

    def get_predecessors(self, v):
        res = []
        for k in self.graph.keys():
            if v in self.graph[k]:
                res.append(k)
        return res

    def get_adjacents(self, v):
        suc = self.get_successors(v)
        pred = self.get_predecessors(v)
        return list(set(suc + pred))

    ## degrees

    def out_degree(self, v):
        return len(self.get_successors(v))

    def in_degree(self, v):
        return len(self.get_predecessors(v))

    def degree(self, v):
        return len(self.get_adjacents(v))

    def all_degrees(self, deg_type = "inout"):
        ''' Computes the degree (of a given type) for all nodes.
        deg_type can be "in", "out", or "inout" '''
        degs = {}
        for v in self.graph.keys():
            match deg_type:
                case "out":
                    degs[v] = self.out_degree(v)
                case "in":
                    degs[v] = self.in_degree(v)
                case "inout":
                    degs[v] = self.degree(v)
        return degs

    def highest_degrees(self, all_deg = None, deg_type = "inout", top = 10):
        if all_deg is None:
            all_deg = self.all_degrees(deg_type)
        ord_deg = sorted(list(all_deg.items()), key=lambda x : x[1], reverse = True)
        return list(map(lambda x:x[0], ord_deg[:top]))


    ## topological metrics over degrees

    def mean_degree(self, deg_type = "inout"):
        ''' average degree of all nodes: sum of all degrees divided by number of nodes'''
        degs = self.all_degrees(deg_type)
        return sum(degs.values())/len(degs)

    def prob_degree(self, deg_type = "inout"):
        # count the number of occurrences of each degree in the network and derive its frequencies
        degs = self.all_degrees(deg_type)
        freqs = {}
        for v in degs.values():
            freqs[v] = freqs.get(v, 0) + 1
        freqs = {k: v / len(degs) for k, v in freqs.items()}
        return freqs

    def print_prob_degree(self, deg_type = "inout"):
        freqs = self.prob_degree(deg_type)
        for k, v in freqs.items():
            print (k, " -> ", v)


    ## BFS and DFS searches

    def reachable_bfs(self, v):
        l = [v]   # list of nodes to be handled
        res = []  # list of nodes to return the result
        while len(l) > 0:
            node = l.pop(0)  # implements a queue: LILO
            if node != v:
                res.append(node)
            for elem in self.graph[node]:
                if elem not in res and elem not in l and elem != node:
                    l.append(elem)
        return res

    def reachable_dfs(self, v):
        l = [v]
        res = []
        while len(l) > 0:
            node = l.pop(0) # implements a stack:
            if node != v:
                res.append(node)
            s = 0
            for elem in self.graph[node]:
                if elem not in res and elem not in l:
                    l.insert(s, elem)
                    s += 1
        return res

    def distance(self, s, d):
        if s == d:
            return 0
        l = [(s,0)]
        visited = [s]
        while len(l) > 0:
            node, dist = l.pop(0)
            for elem in self.graph[node]:
                if elem == d:
                    return dist + 1
                if elem not in visited:
                    l.append((elem, dist + 1))
                    visited.append(elem)
        return None

    def shortest_path(self, s, d):
        if s == d:
            return 0
        l = [(s,[])]
        visited = [s]
        while len(l) > 0:
            node, preds = l.pop(0)
            for elem in self.graph[node]:
                if elem == d:
                    return preds + [node, elem]
                if elem not in visited:
                    l.append((elem, preds + [node]))
                    visited.append(elem)
        return None

    def reachable_with_dist(self, v):
        l = [(v, 0)]
        res = {}
        while len(l) > 0:
            node, dist = l.pop(0)
            if node != v:
                res[node] = dist
            for elem in self.graph[node]:
                if elem not in res and elem not in [x[0] for x in l]:
                    l.append((elem, dist + 1))
        return res

    def mean_distances(self):
        # get all the distances from all nodes to all other nodes
        # and return the mean of all distances
        dists = []
        for v1 in self.graph.keys():
            for v2 in self.graph.keys():
                if v1 != v2:
                    dist = self.distance(v1, v2)
                    if dist:
                        dists.append(dist)
        return sum(dists) / len(dists)


    ## clustering

    def clustering_coef(self, v):
        # get the list of adjancent nodes
        adjs = self.get_adjacents(v)
        if len(adjs) <= 1:
            return 0.0
        # calculate the number of links of the adjacent nodes
        ligs = 0
        # compare pairwisely if nodes in this list are connected between them
        for i in adjs:
            for j in adjs:
                if i != j:
                    # check if i and j are connected to each other; if yes increment link counter
                    if j in self.graph[i]:
                        ligs += 2  # not sure why this should be 2
        return float(ligs) / (len(adjs) * (len(adjs) - 1))

    def all_clustering_coefs(self):
        # go through all the nodes and calculate its cc
        # put those in a dictionary and return
        ccs = {}
        for v in self.graph.keys():
            ccs[v] = self.clustering_coef(v)
        return ccs

    def mean_clustering_coef(self):
        # get all the clustering coefficients
        # and return the mean of all ccs
        ccs = self.all_clustering_coefs()
        return sum(ccs.values()) / len(ccs)


if __name__ == "__main__":
    gr = MyGraph()
    gr.add_vertex(1)
    gr.add_vertex(2)
    gr.add_vertex(3)
    gr.add_vertex(4)
    gr.add_edge(1,2)
    gr.add_edge(2,3)
    gr.add_edge(3,2)
    gr.add_edge(3,4)
    gr.add_edge(4,2)
    gr.print_graph()
    print(gr.size())

    print (gr.get_successors(2))
    print (gr.get_predecessors(2))
    print (gr.get_adjacents(2))

    print (gr.in_degree(2))
    print (gr.out_degree(2))
    print (gr.degree(2))

    print(gr.all_degrees("inout"))
    print(gr.all_degrees("in"))
    print(gr.all_degrees("out"))

    gr2 = MyGraph({1:[2,3,4], 2:[5,6],3:[6,8],4:[8],5:[7],6:[],7:[],8:[]})
    print(gr2.reachable_bfs(1))
    print(gr2.reachable_dfs(1))

    print(gr2.distance(1,7))
    print(gr2.shortest_path(1,7))
    print(gr2.distance(1,8))
    print(gr2.shortest_path(1,8))
    print(gr2.distance(6,1))
    print(gr2.shortest_path(6,1))

    print(gr2.reachable_with_dist(1))

    print(gr.mean_degree())
    print(gr.prob_degree())
    print(gr.mean_distances())
    print (gr.clustering_coef(1))
    print (gr.clustering_coef(2))

    gr3 = MyGraph()
    gr3.graph = {
        'NRAS': ['NF1', 'BRAF'],
        'BRAF': [],
        'NF1': [],
        'ERBB3': [],
        'FLT3': ['NRAS'],
        'PIK3CA': ['FLT3', 'ERBB3', 'PTEN', 'NRAS'],
        'PTEN': [],
        'TP53': ['PIK3CA'],
        'CTNNB1': ['APC', 'PIK3CA'],
        'SMAD4': ['CTNNB1'],
        'NCOR1': ['SMAD4'],
        'APC': [],
        'SF3B1': []
    }
    print(gr3.get_nodes())
    print(gr3.get_edges())
    print(gr3.size())

    print(gr3.out_degree("PIK3CA"))
    print(gr3.in_degree("PIK3CA"))
    print(gr3.degree("PIK3CA"))
    print(gr3.mean_degree(deg_type="inout"))
    print(gr3.mean_degree(deg_type="in"))
    print(gr3.mean_degree(deg_type="out"))
    print(gr3.prob_degree(deg_type="inout"))
    print(gr3.prob_degree(deg_type="in"))
    print(gr3.prob_degree(deg_type="out"))
    gr3.print_prob_degree(deg_type="inout")
    print(gr3.all_clustering_coefs())
    print(gr3.mean_clustering_coef())

    gr4= MyGraph()
    gr4.add_vertex(1)
    gr4.add_vertex(2)
    gr4.add_vertex(3)
    gr4.add_edge(1,2)
    gr4.add_edge(2,3)
    gr4.print_graph()
