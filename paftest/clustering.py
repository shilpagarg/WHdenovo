'''
For dividing the subgraph with weak connections into two
'''
import networkx as nx
import sys

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
    log(len(disconnected), 'disconnected connected components')
    #log(disconnected)
    o = open(out, 'w')
    for chain in disconnected:
        o.write(' '.join(chain) + '\n')
    o.close()
    
    exit()
    
    count = 0
    for sub in disconnected:
        count += 1
        log('start to work on component', count)
        s = nx.Graph(graph.subgraph(sub))
        log('subgraph extracted')
        if canSplit(s):
            log('it can be splitted')
            c2 = 0
            for rem in nx.minimum_edge_cut(s):
                c2 += 1
                log('removing', rem)
                s.remove_edge(rem[0], rem[1])
            log('removed', c2, 'weak edges')
            divided = bfs_connected_components(s)
            for ss in divided:
                fs = nx.Graph(s.subgraph(ss))
                log('outputing it')
                subgraphs.append(fs)
        else:
            log('outputing it')
            subgraphs.append(s)

    return subgraphs
