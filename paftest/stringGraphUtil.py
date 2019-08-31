import itertools
from collections import defaultdict

# Return True if an object is iterable
def isIterable(obj):
    isIterable = True
    try:
        i = iter(obj)
    except TypeError:
        isIterable = False
    return isIterable


# Check that the edge list is sorted
# in ascending order
def edgeListIsSorted(edgeList):
    if len(edgeList) < 2:
        return True
    edgeLens = [edge.len for edge in edgeList]
    isSorted = all(l1 <= l2 for l1,l2 in zip(edgeLens[0:-1], edgeLens[1:]))
    return isSorted

# Reverse an end from 'B' to 'E' or 'E' to 'B'
def reverseEnd(end):
    if end not in ('B', 'E'):
        raise  RuntimeError("end must be one of 'B' or 'E'")
    return 'B' if end=='E' else 'E'


def collectWalksByNodeId(walkList):
    nodeToWalk = defaultdict(list)
    for walk in walkList:
        nodeIds = (n.id for n in walk.nodeIter())
        for nid in nodeIds:
            nodeToWalk[nid].append(walk)
    return nodeToWalk

def printList(l):
    for i in l:
        print(i)

printEdgeList = printList
printNodeList = printList
