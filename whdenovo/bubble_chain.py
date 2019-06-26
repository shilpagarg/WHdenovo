import sys
from collections import defaultdict
import networkx as nx
from itertools import groupby
from vg_reader import vg_read
from util import *

def dfs(graph, start):
	visited, stack = [], [start]
	while stack:
		vertex = stack.pop()
		if vertex not in visited:
			visited.append(vertex)
			stack.extend(graph[vertex] - set(visited))
	return visited

def find_bubble_chains(consec_pairs):
		
	print('the total number of edges %d.' %len(consec_pairs))
	consec_pairs_tmp = defaultdict(set)
	
	for edge in consec_pairs:
		if edge[0] != edge[1]:
			consec_pairs_tmp[edge[0]].add(edge[1])
			consec_pairs_tmp[edge[1]].add(edge[0])

	count = 0
	for i,j in consec_pairs_tmp.items():
		j = list(j)
		if len(j) > 2:
			count += 1
			for k in j:
				if (i, k) in consec_pairs:
					consec_pairs.remove((i, k))
				if (k, i) in consec_pairs:
					consec_pairs.remove((k, i))

	print('the total number of branching points %d.' %count)
	print('the total number of good pairs or edges %d.' %len(consec_pairs))

	# Return all distinct simple paths ending at "targetNode", continuing
	# from "currentPath". "usedNodes" is useful so we can quickly skip
	# nodes we have already added to "currentPath". When a new solution path
	# is found, append it to "answerPaths" and return it.

	nodeToNodess= defaultdict(set)
	nodeToNodes = defaultdict(list)
	for edge in consec_pairs:
		nodeToNodess[edge[0]].add(edge[1])
		nodeToNodess[edge[1]].add(edge[0])

	for i,j in nodeToNodess.items():
		nodeToNodes[i] = list(j)
	start_or_end = set()
	for i,j in nodeToNodes.items():
		if len(nodeToNodes[i]) ==1:
			start_or_end.add(i)
			
	nx_nodeToNodes = nx.MultiDiGraph(nodeToNodes)

	unitigs = []
	hashse = {}
	for i in start_or_end:
		hashse[i]=0

	for i in start_or_end:
		if hashse[i]==0:
			unitigs.append([])
			for j in dfs(nodeToNodess, i):
				unitigs[-1].append(j)
				if j in hashse:
					hashse[j]=1

	print('total number of unitigs', len(unitigs))
	return unitigs

def aux_contigs(unitigs, totalAlnSet):
	unitigs_ends = {}
	unitigs_starts = {}
	unitigs_node_pos = defaultdict(tuple)
	unitigs_enddict = defaultdict(int)
	for var in range(len(unitigs)):
		if len(unitigs[var]) == 1:
			unitigs_ends[unitigs[var][0]] = unitigs[var][0]
			unitigs_starts[unitigs[var][0]] = unitigs[var][0]
			unitigs_enddict[unitigs[var][0]] = 0
		else:
		    unitigs_ends[unitigs[var][0]] = unitigs[var][-1]
		    unitigs_starts[unitigs[var][-1]] = unitigs[var][0]
		    unitigs_enddict[unitigs[var][0]] = 0
		    unitigs_enddict[unitigs[var][-1]] = 0

		for i in range(len(unitigs[var])):
			unitigs_node_pos[unitigs[var][i]] = (var, i)

	aux = aux_unitigs(unitigs)
	for sample in range(len(totalAlnSet.fullReadList)):
		for read in totalAlnSet.fullReadList[sample]:
			unitigDirection = defaultdict(int)
			new_G = []
			new_g = []
			unitig_connection = []
			lastPos = defaultdict(int) # [unitig index, index in unitig], initially it doesn't mean 'last'
			for partial in read.partials:
				for node in partial:
					if node in unitigs_node_pos:
						if len(unitigs[unitigs_node_pos[node][0]]) == 1:
							unitigDirection[unitigs_node_pos[node][0]] = 1
							continue
						currentPos = unitigs_node_pos[node]
						if currentPos[1] > lastPos[currentPos[0]]:
							unitigDirection[currentPos[0]] = 1 # Forward
						elif currentPos[1] < lastPos[currentPos[0]]:
							unitigDirection[currentPos[0]] = -1 # Reverse
						lastPos[currentPos[0]] = currentPos[1]
			for partial in read.partials:
				for node in partial:
					if node in unitigs_ends.keys() or node in unitigs_starts.keys():
						unitigIdx = unitigs_node_pos[node][0]
						direction = unitigDirection[unitigIdx]
						if unitigs_enddict[node] == 0:
							if direction == 0:
								continue
							else:
								# Appending Order Matters!
								if node in unitigs_ends.keys():
									if direction == 1:
										new_g.append(node)
										if unitigs_enddict[unitigs_ends[node]] == 1:
											new_g.append(unitigs_ends[node])
											new_gg = new_g.copy()
											new_G.append(new_gg)
											new_g = []
										else:
											new_g.append(unitigs_ends[node])
										unitigs_enddict[node] = 1
										unitigs_enddict[unitigs_ends[node]] = 1
									elif direction == -1:
										if unitigs_enddict[unitigs_ends[node]] == 1:
											new_gg = new_g.copy()
											new_G.append(new_gg)
											new_g = []
										new_g.append(unitigs_ends[node])
										new_g.append(node)
										unitigs_enddict[unitigs_ends[node]] = 1
										unitigs_enddict[node] = 1
								elif node in unitigs_starts.keys():
									if direction == 1:
										if unitigs_enddict[unitigs_starts[node]] == 1:
											new_gg = new_g.copy()
											new_G.append(new_gg)
											new_g = []
										new_g.append(unitigs_starts[node])
										new_g.append(node)
										unitigs_enddict[unitigs_starts[node]] = 1
										unitigs_enddict[node] = 1
									elif direction == -1:
										new_g.append(node)
										if unitigs_enddict[unitigs_starts[node]] == 1:
											new_g.append(unitigs_starts[node])
											new_gg = new_g.copy()
											new_G.append(new_gg)
											new_g = []
										else:
											new_g.append(unitigs_starts[node])
										unitigs_enddict[unitigs_starts[node]] = 1
										unitigs_enddict[node] = 1
			new_gg = new_g.copy()
			new_G.append(new_gg)
			for new_g in new_G:
				if len(new_g) > 2:
					aux.addConnection(new_g)
				for n in new_g:
					unitigs_enddict[n] = 0

	AUX = aux.returnConfident()
	
	return AUX

def find_contigs(unitigs, AUX):
	consec_pairs_final = defaultdict(set)

	count=0   
	nodeToNodess= defaultdict(set)
	nodeToNodes = defaultdict(list)

	for var in AUX:
		if len(var)>1:
			for i in range(0,len(var)-1):
				if var[i] != var[i+1]:
					nodeToNodess[var[i+1]].add(var[i])
					nodeToNodess[var[i]].add(var[i+1])
	for var in unitigs:
		if len(var) > 1:
			for i in range(0,len(var)-1):
				if var[i] != var[i+1]:
					nodeToNodess[var[i+1]].add(var[i])
					nodeToNodess[var[i]].add(var[i+1])
		if len(var) > 2:
			if var[0] in nodeToNodess[var[len(var)-1]]:
				nodeToNodess[var[len(var)-1]].remove(var[0])
			if var[len(var)-1] in nodeToNodess[var[0]]:
				nodeToNodess[var[0]].remove(var[len(var)-1])
			
	for i,j in nodeToNodess.items():
		nodeToNodes[i] = list(j)
	bubble1 = []
	for key, value in nodeToNodes.items():
		if int(key) < 0:
			count += 1
			bubble1.append(key)

	start_or_end = set()
	for i,j in nodeToNodes.items():
		if len(nodeToNodes[i]) ==1:
			start_or_end.add(i)
	final_ctgs = []
	bubbles = []
	nx_nodeToNodes = nx.MultiDiGraph(nodeToNodes)

	hashse = {}
	for i in start_or_end:
		hashse[i]=0

	for i in start_or_end:
		if hashse[i]==0:
			final_ctgs.append([])
			for j in dfs(nodeToNodess, i):
				final_ctgs[-1].append(j)
				if j in hashse:
					hashse[j]=1
	print("Total number of final blocks",len(final_ctgs))
	fcc = 0
	final_bubble_number = 0
	useful = []
	for c in final_ctgs:
		npc = 0
		for n in c:
			if int(n) < 0:
				npc += 1
				
		if npc >= 2:
			final_bubble_number += npc
			fcc += 1
			useful.append(c)
	print(fcc, 'final blocks have at least 2 bubbles')
	print('with', final_bubble_number, 'bubbles left')
	return useful

def group_readset(contig, totalAlnSet):

	bc_variants = dict()
	phaseg_variant = 0
	for node in contig:
		if int(node) < 0:
			phaseg_variant += 1
			bc_variants[int(node)] = phaseg_variant
	readsets = []
	for sample in range(len(totalAlnSet.bubbleReadMap)):
		readset_ind = set()
		for bubble in list(bc_variants.keys()):
			readset_ind = readset_ind.union(totalAlnSet.bubbleReadMap[sample][bubble])
		readList_ind = []
		for read in readset_ind:
			# Don't partition reads with only one bubbles for now
			readInfo = [read.name]
			for var in read.alleles:
				try:
					readInfo.append((bc_variants[var[0]], var[1]))
				except KeyError:
					pass
			if len(readInfo) < 3:
				continue
			readList_ind.append(readInfo)
		readsets.append(readList_ind)

	return readsets

def bc(locus_file, phase_input_files, t):
	print('Input dataset:', locus_file, phase_input_files)
	totalAlnSet, consec_pairs = vg_read(locus_file, phase_input_files, t)
	unitigs = find_bubble_chains(consec_pairs)
	
	AUX = aux_contigs(unitigs, totalAlnSet)
	final_ctgs = find_contigs(unitigs, AUX)
	total_readsets = []
	for contig in final_ctgs:
		readsets = group_readset(contig, totalAlnSet)
		total_readsets.append(readsets)

	return total_readsets
