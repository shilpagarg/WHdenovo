import sys
import stream
import vg_pb2
from collections import Counter
from collections import defaultdict, OrderedDict
import networkx as nx

unitigs = sys.argv[1]
#normal gam file on bubble space
aux_unitigs=sys.argv[2]

consec_pairs_final = defaultdict(set)
consec_pairs = defaultdict(set)



count=0   
def getAllSimplePaths(originNode, targetNode, nodeToNodes):
	
	return helpGetAllSimplePaths(targetNode,
								 [originNode],
								 set(originNode),
								 nodeToNodes,
								 list(), count)
 
#
# Return all distinct simple paths ending at "targetNode", continuing
# from "currentPath". "usedNodes" is useful so we can quickly skip
# nodes we have already added to "currentPath". When a new solution path
# is found, append it to "answerPaths" and return it.
#
def helpGetAllSimplePaths(targetNode,
						  currentPath,
						  usedNodes,
						  nodeToNodes,
						  answerPaths, count):
	lastNode = currentPath[-1]
	if count > 500:
		count=0
		return []
	if lastNode == targetNode:
		answerPaths.append(list(currentPath))
	else:
		for neighbor in nodeToNodes[lastNode]:
			if neighbor not in usedNodes:
				count+=1
				currentPath.append(neighbor)
				usedNodes.add(neighbor)
				helpGetAllSimplePaths(targetNode,
									  currentPath,
									  usedNodes,
									  nodeToNodes,
									  answerPaths, count)
				usedNodes.remove(neighbor)
				currentPath.pop()
	return answerPaths

nodeToNodess= defaultdict(set)
nodeToNodes = defaultdict(list)

with open(aux_unitigs) as fp:
	for line in fp:
		#g = vg_pb2.Alignment()
		#g.ParseFromString(data)
		var = line.rstrip().split(",")[:-1]
		if len(var)>1:
			for i in range(0,len(var)-1):
				if str(var[i])!=str(var[i+1]):
					nodeToNodess[str(var[i+1])].add(str(var[i]))
					nodeToNodess[str(var[i])].add(str(var[i+1]))

with open(unitigs) as fp:
	for line in fp:
		var = line.rstrip().split(",")[:-1]
		if len(var) > 1:
			for i in range(0,len(var)-1):
				if str(var[i])!=str(var[i+1]):
					nodeToNodess[str(var[i+1])].add(str(var[i]))
					nodeToNodess[str(var[i])].add(str(var[i+1]))
		if len(var) > 2:
			if str(var[0]) in nodeToNodess[str(var[len(var)-1])]:
				nodeToNodess[str(var[len(var)-1])].remove(str(var[0]))
			if str(var[len(var)-1]) in nodeToNodess[str(var[0])]:
				nodeToNodess[str(var[0])].remove(str(var[len(var)-1]))

		
for i,j in nodeToNodess.items():
	if len(j) > 2:
		print (i,j)
	nodeToNodes[i] = list(j)
# count = 0
bubble1 = []
for key, value in nodeToNodes.items():
	if int(key) < 300:
		count += 1
		bubble1.append(key)
print(count)

#print(nodeToNodes.keys())
# keys = []
# with open ('sim.unitigs', 'r') as f:
# 	for line in f:
# 		var = line.rstrip().split(",")
# 		for element in var:
# 			keys.append(element)
# for key in keys:
# 	if key not in nodeToNodes.keys():
# 		print(key)
# 	else:
# 		print('ok')
#print (nodeToNodes)
start_or_end = set()
for i,j in nodeToNodes.items():
	if len(nodeToNodes[i]) ==1:
		start_or_end.add(i)
print('the total number of start or end nodes %d.' %len(start_or_end))
print(start_or_end)
outf = open(sys.argv[3], 'w')
bubbles = []
nx_nodeToNodes = nx.MultiDiGraph(nodeToNodes)
# print all paths between starts and ends
# for i in start_or_end:
# 	for j in start_or_end:
# 		if i!=j and i<j:
# 			if len(getAllSimplePaths(i,j,nodeToNodes))!=0:
# 				#print(getAllSimplePaths(i,j,nodeToNodes))
# 				for k in getAllSimplePaths(i,j,nodeToNodes)[0]:
# 					outf.write(str(k) + ",")
# 					if int(k) < 300:
# 						print(k)
# 						bubbles.append(k)
# 				outf.write("\n")
# 				outf.write("\n")
# 				outf.write("\n")
# 				outf.write("\n")
for i in start_or_end:
	for j in start_or_end:
		if i!=j and int(i)<int(j):
			if nx.has_path(nx_nodeToNodes,i,j) == True:
				print (list(nx.all_simple_paths(nx_nodeToNodes,i,j)))
				if len(list(nx.all_simple_paths(nx_nodeToNodes,i,j)))!=0:
					#print(getAllSimplePaths(i,j,nodeToNodes))
					for k in list(nx.all_simple_paths(nx_nodeToNodes,i,j))[0]:
						outf.write(str(k) + ",")
						if int(k) < 300:
							bubbles.append(k)
					outf.write("\n")
					outf.write("\n")
					outf.write("\n")
					outf.write("\n")
print(bubbles, "\n", len(bubbles))
print(bubble1)
# for element in bubble1:
# 	if element not in bubbles:
# 		print(element)

print("TEST PATHS")
# print(getAllSimplePaths('7370100','7369644', nodeToNodes))
# testgraph = nx.MultiDiGraph(nodeToNodes)
# print('new_Test', nx.shortest_path(testgraph,'7370100','7369644')),
# print('new_test2', nx.shortest_path(testgraph,'7369644','7370100'))
# print('hello')
# print(getAllSimplePaths('7369644','7370100', nodeToNodes))
# print('WHOA')
# path = getAllSimplePaths('7370100','7369644', nodeToNodes)
# test = {}
# if len(path) > 0:
# 	path = path[0]
# 	for i in range(0,len(path)):
# 		test[path[i]] = nodeToNodes[path[i]]
# 		# if len(nodeToNodes[path[i]]) != 2:
# 		# 	print(path[i])
# 		# 	print(nodeToNodes[path[i]])
# 		# if path[i-1] not in nodeToNodes[path[i]] or path[i+1] not in nodeToNodes[path[i]]:
# 		# 	print('gotcha', path[i])
# 		print(path[i], nodeToNodes[path[i]])
# #print('test',getAllSimplePaths('7370100','7369644', nodeToNodes))
# #print('test',getAllSimplePaths('7369644','7370100', nodeToNodes))
# print('test path',getAllSimplePaths('7369644','50', nodeToNodes))
# print('test path',getAllSimplePaths('50','7369644', nodeToNodes))

# #print(test)
# # print(getAllSimplePaths('7370997','7369481', nodeToNodes))
# # print(getAllSimplePaths('7369481','7370997', nodeToNodes))

