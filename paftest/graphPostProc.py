import sys
from collections import defaultdict
from stringGraph import StringGraph, Node, Overlap, NodeCoords, Edge


def log(*content):
    out = []
    for i in content:
        out.append(str(i))
    sys.stderr.write(' '.join(out)+'\n')

gfa = open(sys.argv[1], 'r').read().split('\n')[:-1]
node2node = defaultdict(set)
for line in gfa:
    if line[0] == 'S':
        print(line)
        continue
    if line[0] == 'L':
        tokens = line.split('\t')
        assert len(tokens) == 10
        # format of 'tokens'
        # ['L', 'm54238_180921_173448/43123580/ccs', '+', 'm54328_180924_001027/34734842/ccs', '+', '4149M', '11866', '12684', '0', 'DASH']
        src = tokens[1]
        tgt = tokens[3]
        srcdir = tokens[2]
        tgtdir = tokens[4]
        ovlp = int(tokens[5][:-1])
        srclen = int(tokens[6])
        tgtlen = int(tokens[7])
        nvar = int(tokens[8])
        # 'E' stands for connecting to the original end of 'key' read
        # 'B' stands for connecting to the original beginning of 'key' read
        if srcdir == '+':
            node2node[src].add((ovlp, ovlp/srclen, ovlp/tgtlen, nvar, tgt, 'E', tuple(tokens)))
        else:
            node2node[src].add((ovlp, ovlp/srclen, ovlp/tgtlen, nvar, tgt, 'B', tuple(tokens)))
        if tgtdir == '+':
            node2node[tgt].add((ovlp, ovlp/srclen, ovlp/tgtlen, nvar, src, 'B', tuple(tokens)))
        else:
            node2node[tgt].add((ovlp, ovlp/srclen, ovlp/tgtlen, nvar, src, 'E', tuple(tokens)))

Llinesout = []
seenPair = set()
for src, info in node2node.items():
    if len(info) == 1:
        continue
    info = sorted(list(info), reverse = True)
    #log('selection debugging', src, info)
    maxovlp = info[0][0]
    minovlp = maxovlp * 0.65 # THRESHOLD, for selecting the best ovlpped connections
    good = []
    potential = []
    #bad = []
    dash = []
    for conn in info:
        if conn[3] > 270:
            continue
        if conn[3] > 10:
            good.append(conn)
        elif conn[0] > minovlp and conn[3] > 4 and conn[1] >= 0.55 and conn[2] >= 0.55:
            good.append(conn)
        elif conn[0] > minovlp and conn[3] > 0 and conn[1] >= 0.55 and conn[2] >= 0.55:
            potential.append(conn)
        elif conn[3] == 0 and conn[1] >= 0.4 and conn[2] >= 0.4:
            dash.append(conn)
    if src == 'm64011_181218_235052/56625602/ccs':
        log('for src', src)
        log('good', good)
        log('potential', potential)

    connected = []
    nB, nE, t = 0, 0, 0

    while ( nB < 2 or nE < 2 ) and t < len(good):
        if good[t][5] == 'E':
            if nE <= 2:
                nE += 1
                connected.append(good[t])
        if good[t][5] == 'B':
            if nB <= 2:
                nB += 1
                connected.append(good[t])
        t += 1

    t = 0
    while ( nB < 2 or nE < 2 ) and t < len(potential):
        if potential[t][5] == 'E':
            if nE < 2:
                nE += 1
                connected.append(potential[t])
        if potential[t][5] == 'B':
            if nB < 2:
                nB += 1
                connected.append(potential[t])
        t += 1        

    while ( nB < 2 or nE < 2 ) and t < len(dash):
        if dash[t][5] == 'E':
            if nE < 2:
                nE += 1
                connected.append(dash[t])
        if dash[t][5] == 'B':
            if nB < 2:
                nB += 1
                connected.append(dash[t])
        t += 1   

    #print('N\t'+src+'\tnE:'+str(nE)+'\tnB:'+str(nB))

    #connected.extend(dash)
    #if len(good) == 0:
    #    log('No connection remained for', src)
    for conn in connected:
        
        if (src, conn[4]) not in seenPair and (conn[4], src) not in seenPair:
            #if conn[-2] == 'E':
            Llinesout.append(list(conn[-1]))
            #log(conn[2], src, conn[1], conn[0])
            #elif conn[-2] == 'B':
            #    Llinesout.append(list(conn[-1]))
            #    #log(src, conn[2], conn[1], conn[0])
            #elif conn[3] == 0:
            #    Llinesout.append(list(conn[-1]))
            seenPair.add((src, conn[4]))

def edgeObj2Lline(edge):
    if not isinstance(edge, Edge):
        raise TypeError('Check input object type, I need a stringGraph.Edge object')
    else:
        srcName = edge.src.id
        tgtName = edge.target.id
        if edge.srcEnd == 'E':
            srcDir = '+'
        else:
            srcDir = '-'
        if edge.targetEnd == 'B':
            tgtDir = '+'
        else:
            tgtDir = '-'

        return 'L\t%s\t%s\t%s\t%s\t%s'%(srcName, srcDir, tgtName, tgtDir, str(edge.overlap.getLength())+'M')
def transitive_reduction(Llines):
    
    nodeSet = dict()
    existingPair = set()
    dirdic = {1: True, -1: False}
    G = StringGraph()
    for line in Llines:
        if line[-1] == 'GREEN' or line[-1] == 'ONETWO' or line[-1] == 'DASH':
            srcName = line[1]
            srcLength = int(line[6])
            srcDir = line[2]
            tgtName = line[3]
            tgtLength = int(line[7])
            tgtDir = line[4]
            ovlpL = int(line[5][:-1])
            if srcDir == '+':
                srcStart = srcLength - ovlpL
                srcEnd = srcLength
            else:
                srcStart = 0
                srcEnd = ovlpL
            if tgtDir == '+':
                tgtStart = 0
                tgtEnd = ovlpL
            else:
                tgtStart = tgtLength - ovlpL
                tgtEnd = tgtLength
            #srcName, srcLength, srcStart, srcEnd, srcDir = line[:5]
            #tgtName, tgtLength, tgtStart, tgtEnd, tgtDir = line[5:10]

            if srcDir == '+':
                srcDir = 1
            else:
                srcDir = -1
            if tgtDir == '+':
                tgtDir = 1
            else:
                tgtDir = -1
            
            pair = (srcName, srcDir, tgtName, tgtDir)
            reverse_pair = (tgtName, -1 * tgtDir, srcName, -1 * srcDir)
            if reverse_pair not in existingPair and pair not in existingPair:
                existingPair.add(pair)
                if srcName not in nodeSet:
                    nodeSet[srcName] = Node(srcName, srcLength)
                    G.addNode(nodeSet[srcName])
                if tgtName not in nodeSet:
                    nodeSet[tgtName] = Node(tgtName, tgtLength)
                    G.addNode(nodeSet[tgtName])
                #log('adding Lline', line)
                #log('as ovlp', nodeSet[srcName], srcStart, srcEnd, dirdic[srcDir], nodeSet[tgtName], tgtStart, tgtEnd, dirdic[tgtDir])
                G.addOverlap(Overlap(NodeCoords(nodeSet[srcName], srcStart, srcEnd, dirdic[srcDir]), 
                                     NodeCoords(nodeSet[tgtName], tgtStart, tgtEnd, dirdic[tgtDir])))

    G.transitiveReduce()
    final_L = []
    reducted = set()
    for e in G.edges():
        L = edgeObj2Lline(e)
        l = L.split('\t')
        if (l[1], l[3]) not in reducted and (l[3], l[1]) not in reducted:
            final_L.append(L)
            reducted.add((l[1], l[3]))
            #log(L)

    return final_L

#Llinesout = transitive_reduction(Llinesout)
#log(Llinesout)
for L in Llinesout:
    print('\t'.join(L))
