# Author: Lee Mendelowitz (lmendelo@umiacs.umd.edu)
# Date: 2/20/2013
# A Python implementation of Eugene Myer's String Graph
#
# Reference:
# http://bioinformatics.oxfordjournals.org/content/21/suppl_2/ii79.short 
# DOI: doi: 10.1093/bioinformatics/bti1114
################################################


import copy
import stringGraphUtil as util

####################################################
# Graph Model:

#  - The Nodes "own" their list of edges
#  - The Graph "owns" the nodes
#  - Edges come in pairs, due to bi-directed nature of string graph.
#    Therefore, each edge has a twin with the src/target node swapped.
#  - Edges are computed from an Overlap object, which encodes a sequence
#    overlap between two Node objects.
#  - The Overlap's must be constructed externally. This implementation was built using
#    overlaps computed from Nucmer sequence alignments

#####################
# Define constants
WHITE = 0
RED = 1
BLACK = 2
GREEN = 3

VACANT = WHITE
INPLAY = GREEN
ELIMINATED = RED

FUZZ = 10

####################################################
# Node
class Node(object):

    def __init__(self, id, length):
        self.id = id
        self.length = length
        self.BEdges = []
        self.EEdges = []
        self.CEdges = [] # An edge which implies containment of this Node
        self.NEdges = [] # An edge which implies this Node contains another
        self.BReducedEdges = []
        self.EReducedEdges = []
        self.color = WHITE
        self.edgeIndex = None

    def indexEdges(self):
        self.edgeIndex = {}
        for edge in self.BEdges:
            key = ('B', edge.target.id)
            self.edgeIndex[key] = edge
        for edge in self.EEdges:
            key = ('E', edge.target.id)
            self.edgeIndex[key] = edge

    def addEdge(self, edg):
        if edg.src is not self:
            raise Exception('Adding edge to incorrect node!')

        attr = edg.srcEnd + 'Edges'
        edgeList = self.__dict__[attr]
        edgeList.append(edg)

    def getEdges(self, end):
        assert end in ('B', 'E')
        if end == 'B':
            return self.BEdges
        else:
            return self.EEdges

    # Get an edge by key. The key is a tuple: (end, targetId)
    # Return None if the edge does not exist
    def getEdge(self, key):
        if self.edgeIndex is None:
            self.indexEdges()
        return self.edgeIndex.get(key, None)

    def numOverlapEdges(self):
        return len(self.BEdges) + len(self.EEdges)

    def sortAdjacency(self):
        self.BEdges = sorted(self.BEdges, key = lambda e: e.len)
        self.EEdges = sorted(self.EEdges, key = lambda e: e.len)

    def __str__(self):
        db = len(self.BEdges)
        de = len(self.EEdges)
        outS = '{s.id} len={s.length} deg(B)={db} deg(E)={de}'.format(s=self, db=db, de=de)
        return outS

    # Perform transitive reduction on the set of edges
    # out of one end of this node
    def findTransitiveEdges(self, end):

        reducedEdges = []

        assert(end in 'BE')

        # To match the notation in Myer's algorithm, this node (self) is v
        vwEdges = self.BEdges if end=='B' else self.EEdges

        if not vwEdges:
            return reducedEdges
        #if len(vwEdges) <= 1:
        #    return reducedEdges

        if not util.edgeListIsSorted(vwEdges):
            raise RuntimeError('Adjaceny list must be sorted before transitive reduction')

        # Mark nodes as in play (lines 6-7)
        for vwEdge in vwEdges:
            vwEdge.target.color = INPLAY

        longest = vwEdges[-1].len + FUZZ # The last edge in the list is longest  

        #############################################
        # Lines 9 thru 14
        for vwEdge in vwEdges: # line 9

            w = vwEdge.target
            wEnterEnd = vwEdge.targetEnd # The end we enter w

            if w.color != INPLAY: # line 10
                continue
          
            # If we enter w at 'B', we exit at 'E', and vice versa
            wxEdges = w.EEdges if wEnterEnd == 'B' else w.BEdges
            if not util.edgeListIsSorted(wxEdges):
                raise RuntimeError('Adjaceny list must be sorted before transitive reduction')

            for wxEdge in wxEdges: # line 11

                if (wxEdge.len + vwEdge.len <= longest):
                    x = wxEdge.target
                    if x.color == INPLAY:
                        x.color = ELIMINATED


        ###############################################
        # Lines 15 thru 19:
        for vwEdge in vwEdges: # Line 15
            w = vwEdge.target
            wEnterEnd = vwEdge.targetEnd # The end we enter w

            # If we enter w at 'B', we exit at 'E', and vice versa
            wxEdges = w.EEdges if wEnterEnd == 'B' else w.BEdges
            if not util.edgeListIsSorted(wxEdges):
                raise RuntimeError('Adjaceny list must be sorted before transitive reduction')

            # It's possible that a node X was not marked as ELIMINATED above due to some missing overlap (w_i, x) for some i.
            # This code below tries to handle this case of a missing overlap 
            for edgeRank, wxEdge in enumerate(wxEdges): # line 11
                if ((wxEdge.len < FUZZ) or (edgeRank == 0)):
                   x = wxEdge.target
                   if x.color == INPLAY:
                       # LMM Update: Only mark the edge for reduction if it's placement is within FUZZ of the placement
                       # based on edge (v,x). This is to prevent the algorithm from being overly agressive in marking
                       # an edge as transitive when it probably should not be marked.
                       vxEdge = self.getEdge((end, x.id))
                       assert(vxEdge is not None)
                       vxPlacement = vxEdge.len # placement of the end of x based on edge (v,x)
                       vwxPlacement = vwEdge.len + wxEdge.len # placement of the end of x based on edges (v,w) and (w,x)    
                       # If these placements are consistent, mark the edge (v,x) as transitive.
                       if (abs(vxPlacement - vwxPlacement) < FUZZ):
                           x.color = ELIMINATED # Line 19

        for vwEdge in vwEdges: # Line 20    
            w = vwEdge.target
            if w.color == ELIMINATED:
                reducedEdges.append(vwEdge)
            w.color = VACANT

        return reducedEdges


    ###################################
    # Get the unipath which starts from this node, leaving end
    # Return the list of edges of the path.
    def walkUnipath(self, end):
        assert end in ('B', 'E')

        edges = []
        curNode = self
        curEnd = util.reverseEnd(end)
        startingEnd = curEnd

        while True:
            nextEnd = util.reverseEnd(curEnd)
            edgesOut = curNode.getEdges(nextEnd)

            # Check that this node has one successor
            if len(edgesOut) != 1:
                break

            e = edgesOut[0]

            # Check that the curNode is unique predecessor of the successor
            targetPredecessors = e.target.getEdges(e.targetEnd)
            assert(len(targetPredecessors) > 0)
            if len(targetPredecessors) != 1:
                break

            edges.append(e)
            curNode = e.target
            curEnd = e.targetEnd

            if (curNode.color != GREEN):
                raise RuntimeError('Unipath search is not sane!')

            # If the next node is the starting node, then this is cycle.
            if curNode == self:
                # On a unipath search, the only way to reach the starting node is 
                # to enter it from the same end
                assert(curEnd == startingEnd)
                break

        return edges

    # Find the unipath which includes this node.
    # Return the path as a walk.
    def findUnipath(self):
        forwardEdges = self.walkUnipath('E')
        isCycle = bool(forwardEdges) and (forwardEdges[-1].target == self)

        path = []
        if isCycle:
            # Do not include the last edge of the cycle when reporting the path.
            # This is to avoid duplicating the repeated first/last node when working with the unipath.
            path = forwardEdges[:-1]
        else:
            backwardEdges = self.walkUnipath('B')
            # Reverse the backward edges to orient the path forward through self
            backwardEdges = [e.twin for e in backwardEdges[::-1]]
            path = backwardEdges + forwardEdges

        root = path[0].src if path else self
        walk = Walk(root = root, edges = path)
        walkNodes = walk.nodes()
        return walk


    def reduceTransitiveEdges(self, transitiveEdgeSet):

        irreducibleBEdges = []
        irreducibleEEdges = []

        # Transfer reducible edges to the reduced edge list
        for edge in self.BEdges:
            if edge in transitiveEdgeSet:
                self.BReducedEdges.append(edge)
            else:
                irreducibleBEdges.append(edge)

        for edge in self.EEdges:
            if edge in transitiveEdgeSet:
                self.EReducedEdges.append(edge)
            else:
                irreducibleEEdges.append(edge)

        # Save the irreducible edges
        self.BEdges = irreducibleBEdges
        self.EEdges = irreducibleEEdges


######################################################
# Edge
# A simple wrapper around an overlap object.
# Edges come in pairs (due to the bi-directed nature of string graph)
# The edge is labeled with the sequence of the extension:
# 
# src -------------->
#      target ------------->
#                    ------> label of edge from src to target
# For now, instead of storing the label, we store the length of the overhang
class Edge(object):

    # Make a pair of edges
    @staticmethod
    def makePair(overlap):  
        twin = overlap.makeTwin()
        e1 = Edge(overlap)
        e2 = Edge(twin)
        e1.twin = e2
        e2.twin = e1
        return (e1,e2)

    def __init__(self, overlap):
        self.src = overlap.node1
        self.target = overlap.node2
        self.srcEnd = overlap.node1End # The end edge leaves src
        self.targetEnd = overlap.node2End # The end edge enters target
        self.overlap = copy.copy(overlap) # For now we store a copy of the overlap in the Edge object. This may not be necessary.
        self.twin = None
        self.len = overlap.overhang12
        self.color = WHITE

    # These methods point to the wrapped self.overlap object.
    # These could me made attributes of the edge
    def getSrcOL(self):
        return self.overlap.getLength(1)

    def getTargetOL(self):
        return self.overlap.getLength(2)

    def __str__(self):
        ol = self.overlap.getLength()
        outS = 'SRC:{s.src.id} {s.srcEnd} TGT:[{s.target}] {s.targetEnd} OL:{ol} EXTN:{s.len}'
        return outS.format(s=self, ol=ol)

#####################################################
# Overlap of node1 with node2.
# An overlap consists of two NodeCoords objects, with
# additional attributes computed from the NodeCoords.
class Overlap(object):

    # map from (isPfx, isSfx) pair to overlap code
    endsToCode = { (True, True): 'C',
                   (True, False): 'B',
                   (False, True): 'E',
                   (False, False): 'N' }

    def __init__(self, coords1, coords2):
        self.coords1 = copy.copy(coords1)
        self.coords2 = copy.copy(coords2)
        self.node1 = coords1.node
        self.node2 = coords2.node

        self.resetAttributes()

        # Orient the coordinates to node 1 is forward
        #if (not self.coords1.isForward):
        #    self.coords1.isForward = not self.coords1.isForward
        #    self.coords2.isForward = not self.coords2.isForward

        self.calcAttributes()

    def calcAttributes(self):
        self.computeEnds()
        self.computeOverhang()

    def resetAttributes(self):
        # Attributes set by calcAttributes
        self.node1End = None
        self.node2End = None
        self.type = None # 'BB', 'BE', 'EB', 'EE'

        # ---------------> node 1
        #           --------------> node 2
        # |---------| overhang21
        #                |--------| overhang12
        self.overhang12 = 0
        self.overhang21 = 0

    def getLength(self, id=None):
        if id is None:
            id = 1
        if id == 1:
            return self.coords1.matchedLength()
        else:
            return self.coords2.matchedLength()
        
    # Compute which ends of each read the overlap is
    def computeEnds(self):
        # By YC. W., added two conditions to fit MHC method, when one read is 
        # contained another is prfx/sffx, we still consider this pair.

        # Compute for node 1
        c1 = self.coords1
        self.node1End = Overlap.endsToCode[(c1.isPfx(), c1.isSfx())]
        #print('node1 end type originally is', self.node1End)
        if self.node1End == 'C':
            if c1.isForward:
                #print('C1 is forward, change C to B')
                self.node1End = 'E'
            else:
                #print('C1 is not forward, change C to E')
                self.node1End = 'B'
        #print('node1 end type finally is', self.node1End)

        #print(1, self.node1End)
        # Compute for node 2
        c2 = self.coords2
        self.node2End = Overlap.endsToCode[(c2.isPfx(), c2.isSfx())]
        #print('node2 end type originally is', self.node2End)
        if self.node2End == 'C':
            if c2.isForward:
                #print('C2 is forward, change C to B') 
                self.node2End = 'B'
            else:
                #print('C2 is not forward, change C to E')
                self.node2End = 'E'
        #print(2, self.node2End)
        # Make overlap code ('BB', 'BE', 'EB', or 'EE')
        #print('node2 end type finally is', self.node2End)

        self.type = self.node1End + self.node2End

    # Swap nodes 1 & nodes 2.
    # Orient the overlap so node1 is forward
    def flip(self):
        self.coords2, self.coords1 = self.coords1, self.coords2
        self.node1 = self.coords1.node
        self.node2 = self.coords2.node
        if True:#(not self.coords1.isForward):#(not self.coords1.isForward):  # Changed here. Don't know why it was done like originally
            self.coords1.isForward = not self.coords1.isForward
            self.coords2.isForward = not self.coords2.isForward
        self.calcAttributes()

    def makeTwin(self):
        ovl2 = copy.copy(self)
        ovl2.flip()
        #assert ovl2.coords1.isForward # Changed together with the change above
        return ovl2

    # Compute the overhang of the overlap.
    # Example:
    # -------------> node1
    #        -------------------> node2
    #        |-----| overlap
    #              |------------| overhang
    def computeOverhang(self):
        allowedOLTypes = ('BB', 'BE', 'EB', 'EE')
        if self.type in allowedOLTypes:
            self.overhang21 = self.coords1.unmatchedLength()
            self.overhang12 = self.coords2.unmatchedLength()

#####################################################
# An oriented interval. This could represent the placement of
# a read within a unipath, or the coordinates of a contig overlap.
class Coords(object):
    def __init__(self, start, end, isForward):
        if (end < start):
            raise RuntimeError('Coords: end is less than start!')
        self.start = start
        self.end = end
        self.isForward = isForward

    def length(self):
        return self.end - self.start

    def __str__(self):
        orient = 'Forward' if self.isForward else 'Reverse'
        outS = '({s.start}, {s.end}, {orient})'.format(s = self, orient=orient)
        return outS

#####################################################
# Represents an interval of a node sequence
#
# The indices are given with the python convention.
# start: the starting index, zero based, inclusive
# end: the ending index, exclusive
# If isForward:
#   overlap sequence = nodeSeq[start:end]
# If isReverse:
#   overlap sequence = reverse_complement(nodeSeq[start:end])
#
class NodeCoords(Coords):
    def __init__(self, node, start, end, isForward):
        Coords.__init__(self, start, end, isForward)
        self.node = node

    def isPfx(self):
        return self.start < 20
        return self.start == 0

    def isSfx(self): 
        return self.node.length - self.end < 20
        return self.end == self.node.length

    def isExtreme(self):
        return self.isPfx() or self.isSfx()

    def matchedLength(self):
        return self.length()

    def unmatchedLength(self):
        return self.node.length - self.matchedLength()

    def __str__(self):
        orient = 'Forward' if self.isForward else 'Reverse'
        outS = '({s.node.id}, {s.start}, {s.end}, {orient})'.format(s = self, orient=orient)
        return outS

#####################################################
# Represents the placement of a Node within a walk
#
class NodePlacement(Coords):
    def __init__(self, node, start, end, isForward):
        Coords.__init__(self, start, end, isForward)
        self.node = node

    def __str__(self):  
        orient = 'Forward' if self.isForward else 'Reverse'
        outS = '({s.node.id}, {s.start}, {s.end}, {orient})'.format(s = self, orient=orient)
        return outS

#####################################################
# Represents a walk through a string graph
class Walk(object):

    # keywords: root, edges
    def __init__(self, *args, **kwargs):
        self.root = None
        self.rootIsForward = True # If edges are provided, determine the root oriention from the first edge
        self.edges = []

        if 'root' in kwargs:
            self.root = kwargs['root']
            if not isinstance(self.root, Node):
                raise TypeError('Walk root must be a Node')
            self.edges = kwargs.get('edges', [])
        elif 'edges' in kwargs:
            edges = kwargs['edges']
            if not edges:
                raise RuntimeError('No root provided, and edge list is empty!')
            self.root = edges[0].src
            self.edges = edges

        if not self.checkSane():
            raise RuntimeError('Edge list is not sane!')

        self.computeRootOrientation()

    def checkSane(self):
        if not self.edges:
            return True

        firstEdge = self.edges[0]
        if firstEdge.src != self.root:
            return False

        curNode = self.root
        curEnd = util.reverseEnd(firstEdge.srcEnd)

        for edge in self.edges:
            # Check that this node exits the current node at the appropriate end
            nodeOK = (edge.src == curNode)
            edgeOK = (edge.srcEnd == util.reverseEnd(curEnd))
            if not (edgeOK and nodeOK):
                return False
            curNode = edge.target
            curEnd = edge.targetEnd

        return True

    def computeRootOrientation(self):
        if self.edges:
            self.rootIsForward = (self.edges[0].srcEnd == 'E')
        else:
            self.rootIsForward = True

    def nodeIter(self):
        if not self.edges:
            # We return an iterator over just the single node!
            nbunch = (n for n in (self.root,))
            return (n for n in nbunch)
        srcNodes = (e.src for e in self.edges)
        lastNode = (self.edges[-1].target,)
        return (n for nbunch in (srcNodes, lastNode) for n in nbunch)
        
    def nodes(self):
        return list(self.nodeIter())

    def __str__(self):
        myStr = 'Root: %s\n'%str(self.root)
        myStr += '\n'.join(str(e) for e in self.edges)
        return myStr

    # Compute the layout of the reads in the walk
    # The layout is given as a list of NodePlacements
    def computeLayout(self):

        layout = [NodePlacement(self.root, 0, self.root.length, self.rootIsForward)]
        self.layout = layout

        if not self.edges:
            return layout

        for edge in self.edges:
            curNode = edge.target
            curEnd = layout[-1].end # The last bp position of the last Node in the layout
            newEnd = curEnd + edge.len
            newStart = newEnd - curNode.length
            newStart2 = curEnd - edge.getTargetOL()
            assert(newStart == newStart2)
            isForward = (edge.targetEnd == 'B')
            layout.append(NodePlacement(curNode, newStart, newEnd, isForward))

        self.layout = layout
        return layout

    # Reverse the walk. (Note: this modifies self and does not return a copy)
    def reverse(self):

        if not self.edges:
            self.rootIsForward = not self.rootIsForward
        else:
            self.edges = [e.twin for e in self.edges[::-1]]
            self.root = self.edges[0].src
            self.computeRootOrientation()

#####################################################
class StringGraph(object):
    
    #######################################
    def __init__(self):
        self.nodeMap = {} # Map from node-id to Node

    #######################################
    def addNode(self, node):
        if node.id not in self.nodeMap:
            self.nodeMap[node.id] = node

    #######################################
    def addNodes(self, nodes):
        if util.isIterable(nodes):
            for n in nodes:
                self.addNode(n)
        else:
            self.addNode(n)

    #######################################
    # Add edges
    def addOverlap(self, overlap):
        ovl = overlap
        # Add nodes to graph if necessary
        node1, node2 = ovl.node1, ovl.node2
        self.addNodes((node1, node2))

        # Add the edges to the nodes
        e1, e2 = Edge.makePair(ovl)
        node1.addEdge(e1)
        node2.addEdge(e2)

    #######################################
    def addOverlaps(self, overlaps):
        if util.isIterable(overlaps):
            for ovl in overlaps:
                self.addOverlap(ovl)
        else:
            self.addOverlap(ovl)

    #######################################
    # sort adjacencies for each node in the graph
    def sortAdjacencies(self):
        for nodeid, node in self.nodeMap.items():
            node.sortAdjacency()


    ######################################
    # Use E.W.Myer's transitive reduction algorithm to remove transitive edges
    def transitiveReduce(self):
        self.sortAdjacencies()

        for i, n in self.nodeIter():
            n.color = VACANT

        reducibleEdges = []
        for i, n in self.nodeIter():
            reducibleEdges.extend(n.findTransitiveEdges('B'))
            reducibleEdges.extend(n.findTransitiveEdges('E'))
        reducibleEdgesTwins = [e.twin for e in reducibleEdges]
        reducibleEdges = set(e for el in (reducibleEdges, reducibleEdgesTwins) for e in el)

        for i, n in self.nodeIter():
            n.reduceTransitiveEdges(reducibleEdges)
        
        return reducibleEdges

    def colorNodes(self, color):
        for i, n in self.nodeIter():
            n.color = color

    def colorEdges(self, color):
        for e in self.edgeIter():
            e.color = color

    def findUnipaths(self):
        paths = []
        self.colorNodes(GREEN)
        for i, curNode in self.nodeIter():

            if curNode.color == RED: # Node already placed in path
                continue

            walk = curNode.findUnipath()
            walkNodes = walk.nodes()

            # if everything worked correctly, all the nodes in the path should be GREEN
            # since a node can belong to only one unipath
            walkIsSane = walk.checkSane()
            if not walkIsSane:
                raise RuntimeError('Walk is not sane!')
            nodesOK = all(n.color==GREEN for n in walkNodes)
            if not nodesOK:
                raise RuntimeError('Nodes used in more than one walk!')

            # set all unipath nodes to RED
            for n in walkNodes:
                n.color = RED

            paths.append(walk)

        # Check that all nodes have been assigned to a unipath
        for i, n in self.nodeIter():
            assert n.color == RED

        self.colorNodes(WHITE)

        return paths

    def numNodes(self):
        return len(self.nodeMap)

    def numEdges(self):
        return sum(1 for e in self.edgeIter())

    def nodeIter(self):
        return self.nodeMap.items()

    def edgeIter(self):
        for nid, n in self.nodeIter():
            for edge in n.BEdges:
                yield edge
            for edge in n.EEdges:
                yield edge

    def getNode(self, nodeId):
        return self.__getitem__(nodeId)

    def nodes(self):
        return list(self.nodeIter())

    def edges(self):
        return list(self.edgeIter())

    def checkEdgesAreSane(self):
        for n in self.nodeIter():
            edges = (e for el in (n.BEdges, n.EEdges) for e in el)
            for edge in edges:
                assert(edge.src == n)
                target = edge.target
                targetEnd = edge.targetEnd
                assert(targetEnd in 'BE')
                targetEdges = target.BEdges if targetEnd =='B' else target.EEdges
                if edge.twin not in targetEdges:
                    raise RuntimeError('Graph is not sane! Missing edge: %s'%str(edge.twin))

    # Get a node using the dictionary like syntax:
    # node = StringGraph[nodeId]
    def __getitem__(self, nodeId):
        return self.nodeMap[nodeId]
