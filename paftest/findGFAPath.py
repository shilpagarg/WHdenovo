import argparse

class Graph(object):
    """docstring for Edge"""
    def __init__(self, ):
        self.nodemap = {}
        self.edgeOvlp = {}

    def addEdge(self, node1, node1dir, node2, node2dir, ovlp):
        self.nodemap[node1.name] = node1
        self.nodemap[node2.name] = node2

        if node1dir == '+':
            node1.Enodes.add(node2.name)
        else:
            node1.Bnodes.add(node2.name)
        if node2dir == '+':
            node2.Bnodes.add(node1.name)
        else:
            node2.Enodes.add(node1.name)

        self.edgeOvlp[(node1.name, node2.name)] = ovlp
        self.edgeOvlp[(node2.name, node1.name)] = ovlp

    def getStartOrEndNodes(self):
        self.startOrEnd = set()
        for node, N in self.nodemap.items():
            if N.Bnodes == set() or N.Enodes == set():
                self.startOrEnd.add(node)

    def getRandomLongPath(self, read):
        path = []
        if read != None:
            last = read
        else:
            last = list(self.startOrEnd)[0]
        if self.nodemap[last].Enodes == set():
            now = list(self.nodemap[last].Bnodes)[0]
        else:
            now = list(self.nodemap[last].Enodes)[0]
        path.append(last)
        path.append(now)
        while now not in self.startOrEnd:
            if last in self.nodemap[now].Enodes:
                nx = list(self.nodemap[now].Bnodes)[0]
            else:
                nx = list(self.nodemap[now].Enodes)[0]
            if nx in path:
                print('Cycling!!!!!! Break it at', now)
                break
            path.append(nx)
            last = now
            now = nx
        return path

    def getPathSeqLength(self, path):
        rawLength = 0
        ovlpLength = 0
        for i in range(len(path) - 1):
            rawLength += self.nodemap[path[i]].length
            ovlpLength += self.edgeOvlp[(path[i], path[i+1])]
        rawLength += self.nodemap[path[-1]].length
        seqLength = rawLength - ovlpLength
        return seqLength



class Node(object):
    """docstring for Node"""
    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.Bnodes = set()
        self.Enodes = set()

def readGFA(gfaFIle):
    G = Graph()
    gfa = open(gfaFIle).read().split('\n')[:-1]
    nodesSeq = dict()
    nodes = set()
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
    return G

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('GFA', type = str, metavar = 'GFA', help = 'Input GFA graph file.')
    parser.add_argument('-s', '--start', type = str, metavar = 'READID', help = 'Specify a read id to start a random path.')
    #parser.add_argument('-d', '--depth', type = str, metavar = 'FILE', help = 'Output file from "samtools depth".', required = True)
    #parser.add_argument('-m', '--mincov', type = int, metavar = 'INT', default = 20, help = 'Minimum depth [20] to consider as high-coverage.')
    #parser.add_argument('-min-length', type = int, metavar = 'INT', default = 10000, help = 'Minimum region length [10000] to consider')
    #parser.add_argument('-w', '--window', type = int, default = 5, help = 'While going through each position, allow [5] bp that has coverage lower than threshold.')
    #parser.add_argument('-f', '--format', type = str, metavar = '[BED|region]', default = 'BED', help = 'Output format. BED refers to "chr\\tfrom\\tto" for each line; region refers to "chr:from-to" for each line')
    args = parser.parse_args()
    G = readGFA(args.GFA)
    G.getStartOrEndNodes()
    path = G.getRandomLongPath(args.start)
    ovlpLength = G.getPathSeqLength(path)
    print(path)
    print(ovlpLength)

if __name__ == '__main__':
    main()