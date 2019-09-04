import sys
import networkx as nx
from findGFAPath import Graph, Node

'''
Usage: python this.py GFA listOUT regionLength(bp)
'''
import networkx as nx
import sys

def readGFA(gfaFIle):
    G = Graph()
    gfa = open(gfaFIle).read().split('\n')[:-1]
    nodesSeq = dict()
    nodes = set()
    Llines = []
    for line in gfa:
        if line[0] == 'S':
            tokens = line.split('\t')
            nodesSeq[tokens[1]] = tokens[2]
        elif line[0] == 'L':
            tokens = line.split('\t')
            node1 = tokens[1]
            node1dir = tokens[2]
            node2 = tokens[3]
            node2dir = tokens[4]
            ovlp = int(tokens[5][:-1])
            node1len = int(tokens[6])
            node2lne = int(tokens[7])
            if node1 not in G.nodemap:
                n1 = Node(node1, node1len)
            else:
                n1 = G.nodemap[node1]
            if node2 not in G.nodemap:
                n2 = Node(node2, node1len)
            else:
                n2 = G.nodemap[node2]
            G.addEdge(n1, node1dir, n2, node2dir, ovlp)
            Llines.append([node1, node2])
    return G, Llines

def gfaToNX(Llines):
    g = nx.Graph()
    for line in Llines:
        n1 = line[0]
        n2 = line[1]
        g.add_node(n1)
        g.add_node(n2)
        g.add_edge(n1, n2)
    return g

def bfs_connected_components(graph):
    connected_components = []
    nodes = list(graph.nodes())

    while len(nodes)!=0:
        start_node = nodes.pop()
        queue = [start_node] #FIFO
        visited = [start_node]
        while len(queue)!=0:
            start = queue[0]
            queue.remove(start)
            neighbours = list(graph.neighbors(start))
            #print(neighbours)
            for neighbour in neighbours:
                if neighbour not in visited:
                    queue.append(neighbour)
                    visited.append(neighbour)
                    nodes.remove(neighbour)
        connected_components.append(visited)
        
    return connected_components

def canSplit(graph):
    assert isinstance(graph, nx.Graph)
    for n, c in nx.clustering(graph).items():
        if c <= 0.6:
            return True
    return False
def log(*content):
    out = []
    for i in content:
        out.append(str(i))
    sys.stderr.write(' '.join(out)+'\n')
def divide(graph, out):
    subgraphs = []
    disconnected = bfs_connected_components(graph)
    print(len(disconnected))
    #log('Totally', len(disconnected), 'disconnected connected components')
    o = open(out, 'w')
    for chain in disconnected:
        o.write(' '.join(chain) + '\n')
    o.close()
    starts = [i[-1] for i in disconnected]
    return starts

if len(sys.argv) != 4:
	print('Usage: python this.py GFA listOUT regionLength(bp)')
	exit()
gfa = sys.argv[1]

G, Llines = readGFA(gfa)
nxGraph = gfaToNX(Llines)
starts = divide(nxGraph, sys.argv[2])

regionL = int(sys.argv[3])

G.getStartOrEndNodes()
goodComp = []
for s in starts:
	path = G.getRandomLongPath(s)
	ovlpLength = G.getPathSeqLength(path)
	if ovlpLength > 30000:
		goodComp.append(ovlpLength)
	#print(ovlpLength)
#print(len(goodComp))
#print('%d long components, with total length %d' % (len(goodComp), sum(goodComp)), sum(goodComp)/regionL, 'times of original length')
