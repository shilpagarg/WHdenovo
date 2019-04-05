import sys
import stream
import vg_pb2
from collections import Counter
from collections import defaultdict, OrderedDict
import networkx as nx

locus_file = sys.argv[1]
gam_file=sys.argv[2]
consec_pairs = defaultdict(set)

def vg_reader(locus_file, gam_file):
	"""
	input: sorted locus and sorted GAM file output from vg.
	output: sorted readset for core DP.
	assumptions: 
	1. locus file consists of linear ordering of simple bubbles only and hence sorted. Each locus file does not contain start and end vertex.
	2. paths in the locus should be covered by atleast one pacbio read.
	2. GAM file is sorted and restricted to locus file.
	3. files consists of all DAG connected components.
	4. add variant only when it identifies the branch uniquely.
	"""

	locus_count=0
	prev_startsnarl = 0
	prev_endsnarl = 0
	locus_branch_mapping=OrderedDict()
	locus_count=0
	prev_startsnarl = 0
	prev_startsnarl_orientation = -1
	prev_endsnarl = 0
	prev_endsnarl_orientation = -1
	reads_dict = defaultdict(list)
	with stream.open(str(locus_file), "rb") as istream:
		for data in istream:
			l = vg_pb2.SnarlTraversal()
			l.ParseFromString(data)
			#TODO: make ordered doctionary locus_branch_mapping
			# handle forward and backward case of nodes
			current_startsnarl = l.snarl.start.node_id
			current_startsnarl_orientation = l.snarl.start.backward
			current_endsnarl = l.snarl.end.node_id
			current_endsnarl_orientation = l.snarl.end.backward
			path_in_bubble =[]

			if len(l.visits) ==0:
				#TODO: for now, assumed, all nodes in path are either forward or backward
				if l.snarl.start.backward == True:
					path_in_bubble.append(tuple ((l.snarl.end.node_id,l.snarl.start.node_id)))
				else:
					path_in_bubble.append(tuple ((l.snarl.start.node_id,l.snarl.end.node_id)))
			else:
				#TODO: for now, assumed, all nodes in path are either forward or backward
				if l.snarl.start.backward == True:
					path_in_bubble.append(tuple ((l.snarl.end.node_id, l.visits[-1].node_id)))
					for i in range(0,len(l.visits)-1):
						path_in_bubble.append(tuple((l.visits[i+1].node_id, l.visits[i].node_id)))
					path_in_bubble.append(tuple ((l.visits[0].node_id,l.snarl.start.node_id)))
				else:
					path_in_bubble.append(tuple ((l.snarl.start.node_id,l.visits[0].node_id)))
					for i in range(0,len(l.visits)-1):
						path_in_bubble.append(tuple((l.visits[i].node_id, l.visits[i+1].node_id)))
					path_in_bubble.append(tuple ((l.visits[-1].node_id, l.snarl.end.node_id))) 

			if current_startsnarl == prev_startsnarl and current_endsnarl == prev_endsnarl and current_endsnarl_orientation == prev_endsnarl_orientation and prev_startsnarl_orientation == current_startsnarl_orientation:
				per_locus.append(path_in_bubble)
			else:
				locus_count=locus_count-1
				per_locus = []
				per_locus.append(path_in_bubble)
			prev_startsnarl = current_startsnarl
			prev_startsnarl_orientation = current_startsnarl_orientation
			prev_endsnarl = current_endsnarl
			prev_endsnarl_orientation = current_endsnarl_orientation
			locus_branch_mapping[locus_count]=per_locus

	#print("LOVCUS BRANCH MAPPING", locus_branch_mapping)
	print('The number of hets:')
	het_count= 0
	for k,v in locus_branch_mapping.items():
		if len(v) > 1:
			het_count = het_count +1
	print(het_count)
	# keep branch of paths in each bubble.
	alleles_per_pos= defaultdict()
	for k,v in locus_branch_mapping.items():
		alleles_per_pos[k]=len(v)

	# both simple and complex bubbles: key is the values in locus_branch_mapping and value is triplet(locus, branch, alleles)
	reverse_mapping= defaultdict(list)
	for k,v in locus_branch_mapping.items():
		if len(v) > 1: # more than one branch
			for i,b in enumerate(v):
				if len(b) > 0:
					for p,j in enumerate(b):
						reverse_mapping[j].append([k,i, len(v)]) # in complex bubbles, a node can map to multiple branches.


	# both simple and complex bubbles: extract reads from GAM file associated with the locus and create a sorted readset.
	# in complex bubble, set of nodes uniquely determine the path. 
	count =0
	duplicated = 0
	#TODO: consider reads with only positive score.
	with stream.open(str(gam_file), "rb") as istream:
		for data in istream:
			g = vg_pb2.Alignment()
			g.ParseFromString(data) 
			# hard-coded source id, mapping quality and other values.
			val1 = True
			val2 = False

			count1 =0
			count2=0
			score = g.score/len(g.sequence)

			#if score > 0.2:
			#	continue
			read= [] # create read for each read alignment
			prev_tmp=[]
			prev_locus= -1
			locus = -1
			#print(g.name)
			for i in range(0,len(g.path.mapping)-1):
			#for i in g.path.mapping: # go over the mapping in a read
				# TODO: check for forward or reverse strand, we may not need it for DAG.
				edge1 = tuple((int(g.path.mapping[i].position.node_id), int(g.path.mapping[i+1].position.node_id))) # go over nodes in a mapping
				edge2 = tuple((int(g.path.mapping[i+1].position.node_id), int(g.path.mapping[i].position.node_id))) # go over nodes in a mapping
#				print(edge1, edge2)
				if edge1 in reverse_mapping or edge2 in reverse_mapping: # handle start and sink node.
					if edge1 in reverse_mapping:
#						print('in reverse mapping')
						qualities = [10]* reverse_mapping[edge1][0][2]
						node_inf= [tuple(i[0:2]) for i in reverse_mapping[edge1]] # consider (locus, branch)
					else:
#						print('in reverse mapping')
						qualities = [10]* reverse_mapping[edge2][0][2]
						node_inf= [tuple(i[0:2]) for i in reverse_mapping[edge2]]
					tmp = [x for x in node_inf]
					if prev_locus != tmp[0][0]:
						prev_tmp = tmp
						prev_locus = tmp[0][0]
						
					interset_tmp= list(set(tmp).intersection(set(prev_tmp)))
					if len(prev_tmp) > 0 and len(set(tmp).intersection(set(prev_tmp)))==1: # for complicated bubbles, but with Top-k paths. combination of some nodes uniquely determine branch.
						qualities[interset_tmp[0][1]] = 0
						if i== len(g.path.mapping)-2:
							#read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)
							reads_dict[g.name+"_"+str(g.query_position)].append(interset_tmp[0][0])
							read.append(interset_tmp[0][0])
							#print('i am in bubble')
						else:
							next_edge1 = tuple((int(g.path.mapping[i+1].position.node_id), int(g.path.mapping[i+2].position.node_id)))
							next_edge2 = tuple((int(g.path.mapping[i+2].position.node_id), int(g.path.mapping[i+1].position.node_id)))

							if next_edge1 not in reverse_mapping and next_edge2 not in reverse_mapping:
								#read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)	
								reads_dict[g.name+"_"+str(g.query_position)].append(interset_tmp[0][0])
								read.append(interset_tmp[0][0])
								#print('i am in bubble')
						locus= interset_tmp[0][0]
				else:
					read.append(int(g.path.mapping[i].position.node_id))
					read.append(int(g.path.mapping[i+1].position.node_id))
					reads_dict[g.name+"_"+str(g.query_position)].append(int(g.path.mapping[i].position.node_id))
					reads_dict[g.name+"_"+str(g.query_position)].append(int(g.path.mapping[i+1].position.node_id))
			#print(read)

			# for every pair of bubbles or bubble-node
			for k in range(0,len(read)-1):
				pair1= str(read[k])+"_"+str(read[k+1]) # not taking care of reverse direction now
				pair2= str(read[k+1])+"_"+str(read[k])
				# should take of direction, not adding pairs reverse of each other
				if pair2 in consec_pairs:
					consec_pairs[pair2].add(g.name)
				else:
					consec_pairs[pair1].add(g.name)

	return reads_dict, consec_pairs

reads_dict, cosec_pairs = vg_reader(locus_file, gam_file)
print('i am here')
print('the total number of edges %d.' %len(consec_pairs))
consec_pairs_final=defaultdict(set)

# pairs to trust
for i,j in consec_pairs.items():
	#if len(j) > 0 and len(j) < 45:
	# if len(j) > 8 and len(j) < 200:
	if len(j) > 5:
		consec_pairs_final[i] = j

print('the total number of edges %d after removing errorenous ones.' %len(consec_pairs_final))
#print(consec_pairs_final)
consec_pairs_tmp=defaultdict(set)
consec_pairs_tmp_end=defaultdict(set)
start_or_end = set()
# find the points of branches by also considering direction of reads
for i,j in consec_pairs_final.items():
	start = int(i.split("_")[0])
	end = int(i.split("_")[1])
	if start!=end:
		consec_pairs_tmp[start].add(end)
		consec_pairs_tmp[end].add(start)

#print(consec_pairs_tmp)
count = 0
for i,j in consec_pairs_tmp.items():
	if len(j) > 2:
		count+=1
		print(i,j)
		for k in j:
			if str(i)+"_"+str(k) in consec_pairs_final:
				del consec_pairs_final[str(i)+"_"+str(k)]
			if str(k)+"_"+str(i) in consec_pairs_final:
				del consec_pairs_final[str(k)+"_"+str(i)]

#for i,j in consec_pairs_tmp_end.items():
#	if len(j) >= 2:
#		count+=1
#		for k in j:
#			if str(i)+"_"+str(k) in consec_pairs_final:
#				del consec_pairs_final[str(k)+"_"+str(i)]
print('the total number of branching points %d.' %count)
print('the total number of good pairs or edges %d.' %len(consec_pairs_final))

#print(consec_pairs_final)

from collections import defaultdict 
   
def getAllSimplePaths(originNode, targetNode, nodeToNodes):
	return helpGetAllSimplePaths(targetNode,
								 [originNode],
								 set(originNode),
								 nodeToNodes,
								 list())
 
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
						  answerPaths):
	lastNode = currentPath[-1]
	if lastNode == targetNode:
		answerPaths.append(list(currentPath))
	else:
		for neighbor in nodeToNodes[lastNode]:
			if neighbor not in usedNodes:
				currentPath.append(neighbor)
				usedNodes.add(neighbor)
				helpGetAllSimplePaths(targetNode,
									  currentPath,
									  usedNodes,
									  nodeToNodes,
									  answerPaths)
				usedNodes.remove(neighbor)
				currentPath.pop()
	return answerPaths

nodeToNodess= defaultdict(set)
nodeToNodes = defaultdict(list)
for i,j in consec_pairs_final.items():
	if i.split("_")[0] != i.split("_")[1]:
		nodeToNodess[i.split("_")[0]].add(i.split("_")[1])
		nodeToNodess[i.split("_")[1]].add(i.split("_")[0])

for i,j in nodeToNodess.items():
	nodeToNodes[i] = list(j)

#print(nodeToNodes)
for i,j in nodeToNodes.items():
	if len(nodeToNodes[i]) ==1:
		start_or_end.add(i)
		
print('the total number of start or end nodes %d.' %len(start_or_end))
print(start_or_end)
nx_nodeToNodes = nx.MultiDiGraph(nodeToNodes)
outf = open(sys.argv[3], 'w')
# print all paths between starts and ends
# for i in start_or_end:
# 	for j in start_or_end:
# 		if i!=j and int(i) < int(j):
# 			if len(getAllSimplePaths(i,j,nodeToNodes))!=0:
# 				for k in getAllSimplePaths(i,j,nodeToNodes)[0]:
# 					outf.write(str(k) + ",")
# 				outf.write("\n")

def dfs(graph, start):
	visited, stack = [], [start]
	while stack:
		vertex = stack.pop()
		if vertex not in visited:
			visited.append(vertex)
			stack.extend(graph[vertex] - set(visited))
	return visited

hashse = {}
for i in start_or_end:
	hashse[i]=0

for i in start_or_end:
	if hashse[i]==0:
		for j in dfs(nodeToNodess, str(i)):
			outf.write(str(j) + ",")
			if j in hashse:
				hashse[j]=1
		outf.write("\n")
	

#for i in start_or_end:
#	for j in start_or_end:
#		if i!=j and int(i) < int(j):
#			if nx.has_path(nx_nodeToNodes,i,j) == True:
				#print (nx.shortest_path(nx_nodeToNodes,i,j))
#				if len(list(nx.all_simple_paths(nx_nodeToNodes,i,j)))!=0:
					#print(getAllSimplePaths(i,j,nodeToNodes))
#					for k in list(nx.all_simple_paths(nx_nodeToNodes,i,j))[0]:
#						outf.write(str(k) + ",")
#					outf.write("\n")

