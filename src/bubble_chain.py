import sys
import stream
import vg_pb2
from collections import Counter
from collections import defaultdict, OrderedDict
import networkx as nx
from itertools import groupby
#from vg_reader_simplified import vg_reader2
import linecache
from vg_reader import vg_reader_multithreading
import time
from multiprocessing import Pool, Queue, Manager
from ctypes import c_char_p
from util import *

def dfs(graph, start):
	visited, stack = [], [start]
	while stack:
		vertex = stack.pop()
		if vertex not in visited:
			visited.append(vertex)
			stack.extend(graph[vertex] - set(visited))
	return visited

def nodeCoverage(node, dictionary):
	if node < 0:
		c = 0
		for ind in dictionary:
			c += len(ind[node])
		try:
			c_child = len(dictionary[2][node])
		except IndexError:
			c_child = len(dictionary[0][node])
		return c, c_child
	else:
		c = dictionary[node]
		return c

def nodeCoverage_hard(node, dictionary):
        if node < 0:
                c = 0
                for node in dictionary:
                        c += len(dictionary[node])
        else:
                c = len(dictionary[node])
        return c

def find_bubble_chains(consec_pairs, nodeCoveredByRead, bubbleCoverByRead, bubble_of_interest):
		
	print('the total number of edges %d.' %len(consec_pairs))
	t = 0
	#for k, v in bubbleCoverByRead[2].items():
	#	if len(v) > 80:
	#		t += 1
	#	print(k, len(v))
	consec_pairs_final = set()
	#consec_pairs_final=defaultdict(int)
	#consec_pairs_final=defaultdict(set)
	# pairs to trust

	node2nodes = defaultdict(set)
	for i in consec_pairs:
		if i.split("_")[0] != i.split("_")[1]:
			node2nodes[int(i.split("_")[0])].add(int(i.split("_")[1]))
			node2nodes[int(i.split("_")[1])].add(int(i.split("_")[0]))



	for i in consec_pairs:
	#for i,j in consec_pairs.items():
		#if len(j) > 0 and len(j) < 45:
		# if len(j) > 8 and len(j) < 200:
		node1 = int(i.split('_')[0])
		node2 = int(i.split('_')[1])
		c_child1 = None
		c_child2 = None
		if node1 < 0:
			c1, c_child1 = nodeCoverage(node1, bubbleCoverByRead)
		else:
			c1 = nodeCoverage(node1, nodeCoveredByRead)
		if node2 < 0:
			c2, c_child2 = nodeCoverage(node2, bubbleCoverByRead)
		else:
			c2 = nodeCoverage(node2, nodeCoveredByRead)
		
		if c1 > 3000000 or c1 < 1:
			try:
				del node2nodes[node1]
			except:
				pass
			if node2 in node2nodes.keys():
				node2nodes[node2].remove(node1)
		if c2 > 3000000 or c2 < 1:
			try:
				del node2nodes[node2]
			except:
				pass
			if node1 in node2nodes.keys():
				node2nodes[node1].remove(node2)

		if c1 > 3000000 or c1 < 1 or c2 > 3000000 or c2 < 1:
			#print('removing edge',node1, c1, node2, c2, 'child coverage', c_child1, c_child2)
			continue
		else:
			consec_pairs_final.add(i)
			#consec_pairs_final[i] = j
			#print('keeping edge', node1, c1, node2, c2, 'child coverage', c_child1, c_child2)


		'''
		if len(j) >= 5 and len(j) <= 120:
			consec_pairs_final[i] = j
		elif len(j) < 5:
			node1 = int(i.split('_')[0])
			node2 = int(i.split('_')[1])
			node1good, node2good = False, False
			if node1 < 0:
				c1 = nodeCoverage(node1, bubbleCoverByRead)
				if c1 >= 5 and c1 <= 120:
					node1good = True
			else:
				c1 = nodeCoverage(node1, nodeCoveredByRead)
				if c1 >= 5 and c1 <= 120:
					node1good = True
			if node2 < 0:
				c2 = nodeCoverage(node2, bubbleCoverByRead)
				if c2 >= 5 and c2 <= 120:
                                        node2good = True
			else:
				c2 = nodeCoverage(node2, nodeCoveredByRead)
				if c2 >= 5 and c2 <= 120:
					node2good = True
			print('bad coverage, checking each side', c1, c2)
			if node1good and node2good:
				consec_pairs_final[i] = j
		'''

	print('the total number of edges %d after removing errorenous ones.' %len(consec_pairs_final))
	consec_pairs_tmp = defaultdict(set)
	start_or_end = set()
	# find the points of branches by also considering direction of reads
	for i in consec_pairs_final:
	#for i,j in consec_pairs_final.items():
		start = int(i.split("_")[0])
		end = int(i.split("_")[1])
		if start!=end:
			consec_pairs_tmp[start].add(end)
			consec_pairs_tmp[end].add(start)

	count = 0
	removed = defaultdict(set)
	for i,j in consec_pairs_tmp.items():
		j = list(j)
		if len(j) > 2:
			count += 1
			j = sorted(j)
			for k in j[1:]:  # instead of for k in j, randomly maintain one edge here.
				if str(i)+"_"+str(k) in consec_pairs_final:
					consec_pairs_final.remove(str(i)+"_"+str(k))
					#del consec_pairs_final[str(i)+"_"+str(k)]
				if str(k)+"_"+str(i) in consec_pairs_final:
					consec_pairs_final.remove(str(k)+"_"+str(i))
					#del consec_pairs_final[str(k)+"_"+str(i)]
				#print('removing branching edge:',i,k)

	print('the total number of branching points %d.' %count)
	print('the total number of good pairs or edges %d.' %len(consec_pairs_final))

	# Return all distinct simple paths ending at "targetNode", continuing
	# from "currentPath". "usedNodes" is useful so we can quickly skip
	# nodes we have already added to "currentPath". When a new solution path
	# is found, append it to "answerPaths" and return it.
	#	

	nodeToNodess= defaultdict(set)
	nodeToNodes = defaultdict(list)
	for i in consec_pairs_final:
	#for i,j in consec_pairs_final.items():
		if i.split("_")[0] != i.split("_")[1]:
			nodeToNodess[int(i.split("_")[0])].add(int(i.split("_")[1]))
			nodeToNodess[int(i.split("_")[1])].add(int(i.split("_")[0]))

	for i,j in nodeToNodess.items():
		nodeToNodes[i] = list(j)

	for i,j in nodeToNodes.items():
		if len(nodeToNodes[i]) ==1:
			start_or_end.add(i)
			
	#print('the total number of start or end nodes %d.' %len(start_or_end))
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
	for k, v in node2nodes.items():
		if len(v) == 0:
			unitigs.append([k])
	print('total number of unitigs', len(unitigs))
	for i in unitigs:
		print(i)
	return unitigs

def aux_contigs(unitigs, reads_dict, bubble_of_interest):
	unitigs_ends = {}
	unitigs_starts = {}
	#bubble_count = 0
	#bubble_dict = {}
	unitigs_node_pos = defaultdict(tuple)
	single = 0
	for var in range(len(unitigs)):
		if len(unitigs[var]) == 1:
			if unitigs[var][0] < 0:
				single += 1
			unitigs_ends[unitigs[var][0]] = unitigs[var][0]
			unitigs_starts[unitigs[var][0]] = unitigs[var][0]
		else:
		    unitigs_ends[unitigs[var][0]] = unitigs[var][-1]
		    unitigs_starts[unitigs[var][-1]] = unitigs[var][0]
		for i in range(len(unitigs[var])):
			unitigs_node_pos[unitigs[var][i]] = (var, i)
		#for i in range(len(var)):
		#	if int(var[i]) < 0:
		#		bubble_count += 1
		#		if var[i] in bubble_dict.keys():
		#			bubble_dict[var[i]] += 1
		#		else:
		#			bubble_dict[var[i]] = 1
	#print('DEBUG::unitigs_node_pos', unitigs_node_pos)
	unitigs_enddict = defaultdict(int)

	for i in unitigs_ends.keys():
		unitigs_enddict[i] = 0
	for i in unitigs_starts.keys():
		unitigs_enddict[i] = 0

	#nodes_list = set()
	#dummy_list = ['0']*10
	orderalignment = defaultdict(list)
	orderalignment = defaultdict(lambda: [-1]*12, orderalignment)

	for i,j in reads_dict.items():
		canu_name = '_'.join(i.split("_")[:-1])
		canu_chunk_num = int(i.split("_")[-1])
		orderalignment[canu_name].insert(canu_chunk_num, j)

	new_orderalignment = defaultdict(list)
	for k,v in orderalignment.items():
		new_orderalignment[k] = [x for x in v if x != -1]

	norder = defaultdict(list)
	for k,v in new_orderalignment.items():
		norder[k] = [x[0] for x in groupby(v)]

	new_orderalignment = norder

	#aux_unitigs = []
	#bubble_dict2 = dict.fromkeys(bubble_dict.keys(), 0)
	inside_aux = set()
	connected_pair = set()
	#print('unitig_ends', unitigs_ends)
	#print('unitig_starts', unitigs_starts)
	aux = aux_unitigs(unitigs)
	for k,v in new_orderalignment.items():
		new_G = []
		new_g = []
		unitig_connection = []
		unitigDirection = defaultdict(int)
		lastPos = defaultdict(int) # [unitig index, index in unitig], initially it doesn't mean 'last'
		for i in range(0,len(v)):
			for j in range(0,len(v[i])):
				if v[i][j] in unitigs_node_pos:
					if len(unitigs[unitigs_node_pos[v[i][j]][0]]) == 1:
						unitigDirection[unitigs_node_pos[v[i][j]][0]] = 1
						continue
					currentPos = unitigs_node_pos[v[i][j]]
					if currentPos[1] > lastPos[currentPos[0]]:
						unitigDirection[currentPos[0]] = 1 # Forward
					elif currentPos[1] < lastPos[currentPos[0]]:
						unitigDirection[currentPos[0]] = -1 # Reverse
					lastPos[currentPos[0]] = currentPos[1]
		#print('\nLooking at global alignment of', k, v)
		#print('\nunitigs_enddict begin', unitigs_enddict)
		#print('unitigDirection', unitigDirection)
		for i in range(0,len(v)):
			for j in range(0,len(v[i])):
				if v[i][j] in unitigs_ends.keys() or v[i][j] in unitigs_starts.keys():
					unitigIdx = unitigs_node_pos[v[i][j]][0]
					direction = unitigDirection[unitigIdx]
					#print('start to look at node:',v[i][j])
					if unitigs_enddict[v[i][j]] == 0:
						if direction == 0:
							continue
						else:
							# Appending Order Matters!
							if v[i][j] in unitigs_ends.keys():
								if direction == 1:
									new_g.append(v[i][j])
									if unitigs_enddict[unitigs_ends[v[i][j]]] == 1:
										new_g.append(unitigs_ends[v[i][j]])
										new_gg = new_g.copy()
										new_G.append(new_gg)
										new_g = []
									else:
										new_g.append(unitigs_ends[v[i][j]])
									unitigs_enddict[v[i][j]] = 1
									unitigs_enddict[unitigs_ends[v[i][j]]] = 1
								elif direction == -1:
									if unitigs_enddict[unitigs_ends[v[i][j]]] == 1:
										new_gg = new_g.copy()
										new_G.append(new_gg)
										new_g = []
									new_g.append(unitigs_ends[v[i][j]])
									new_g.append(v[i][j])
									unitigs_enddict[unitigs_ends[v[i][j]]] = 1
									unitigs_enddict[v[i][j]] = 1
							elif v[i][j] in unitigs_starts.keys():
								if direction == 1:
									if unitigs_enddict[unitigs_starts[v[i][j]]] == 1:
										new_gg = new_g.copy()
										new_G.append(new_gg)
										new_g = []
									new_g.append(unitigs_starts[v[i][j]])
									new_g.append(v[i][j])
									unitigs_enddict[unitigs_starts[v[i][j]]] = 1
									unitigs_enddict[v[i][j]] = 1
								elif direction == -1:
									new_g.append(v[i][j])
									if unitigs_enddict[unitigs_starts[v[i][j]]] == 1:
										new_g.append(unitigs_starts[v[i][j]])
										new_gg = new_g.copy()
										new_G.append(new_gg)
										new_g = []
									else:
										new_g.append(unitigs_starts[v[i][j]])
									unitigs_enddict[unitigs_starts[v[i][j]]] = 1
									unitigs_enddict[v[i][j]] = 1
		new_gg = new_g.copy()
		new_G.append(new_gg)
		#print('NewG here:', new_G)
		#print('unitigs_enddict iterhou', unitigs_enddict)
		for new_g in new_G:
			if len(new_g) > 2:
				aux.addConnection(new_g)
				#if new_g[0] not in inside_aux or new_g[-1] not in inside_aux:
				#if new_g[0] not in inside_aux:
				#		unitigs_enddict[new_g[0]] = 0
					#if new_g[-1] not in inside_aux:
				#		unitigs_enddict[new_g[-1]] = 0
					#if bubble_of_interest in new_g:
					#	print(k, v, 'aux adding', new_g)
				#aux_unitigs.append(new_g)
					#print('Adopted!')
				#for n in new_g[1:-1]:
				#	inside_aux.add(n)
			#else:
			for n in new_g:
				#if n not in inside_aux:
				unitigs_enddict[n] = 0
					#try:
					#	unitigs_enddict[unitigs_ends[n]] = 0
					#except KeyError:
					#	unitigs_enddict[unitigs_starts[n]] = 0
		#print('unitigs_enddict END', unitigs_enddict)
	#print('DUBUG::aux_unitigs', aux_unitigs)
	AUX = aux.returnConfident()
	print('Final AUX', AUX)
	return AUX

def find_contigs(unitigs, AUX):
	consec_pairs_final = defaultdict(set)

	count=0   
	nodeToNodess= defaultdict(set)
	nodeToNodes = defaultdict(list)

	for var in AUX:
		if len(var)>1:
			for i in range(0,len(var)-1):
				if str(var[i])!=str(var[i+1]):
					nodeToNodess[str(var[i+1])].add(str(var[i]))
					nodeToNodess[str(var[i])].add(str(var[i+1]))
	for var in unitigs:
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
	#print('the total number of start or end nodes %d.' %len(start_or_end))
	final_ctgs = []
	bubbles = []
	nx_nodeToNodes = nx.MultiDiGraph(nodeToNodes)

	hashse = {}
	for i in start_or_end:
		hashse[i]=0

	for i in start_or_end:
		if hashse[i]==0:
			final_ctgs.append([])
			for j in dfs(nodeToNodess, str(i)):
				final_ctgs[-1].append(j)
				if j in hashse:
					hashse[j]=1
	print("Total number of final blocks",len(final_ctgs))
	fcc = 0
	final_bubble_number = 0
	for c in final_ctgs:
		has_bubble = False
		npc = 0
		print(c)
		for n in c:
			if int(n) < 0:
				npc += 1
				final_bubble_number += 1
		if npc >= 2:
			fcc += 1
			#print(c)
	print(fcc, 'usable final blocks')
	print('with', final_bubble_number, 'bubbles left')
	return final_ctgs



def group_readest(contig, bubbleCoverByRead, readnames, readdict_allele):

	bc_variants = dict()
	phaseg_variant = 0
	for node in contig:
		if int(node) < 0:
			phaseg_variant += 1
			bc_variants[int(node)] = phaseg_variant
	
	readsets = []
	for ind in range(len(readdict_allele)):

		readset_ind_idx = set()
		for b in list(bc_variants.keys()):
			readset_ind_idx = readset_ind_idx.union(bubbleCoverByRead[ind][b])
		readnameset = set()
		for i in readset_ind_idx:
			readnameset.add(readnames[ind][i-1])

		orderalignment_allele = defaultdict(list)
		orderalignment_allele = defaultdict(lambda: [-1]*12, orderalignment_allele)

		for i,j in readdict_allele[ind].items():
			canu_name = '_'.join(i.split("_")[:-1])
			if canu_name in readnameset:
				canu_chunk_num = int(i.split("_")[-1])
				orderalignment_allele[canu_name].insert(canu_chunk_num, j)

		new_orderalignment_allele = defaultdict(list)
		for k,v in orderalignment_allele.items():
			new_orderalignment_allele[k] = [x for x in v if x != -1]

		norder_allele = defaultdict(list)
		for k,v in new_orderalignment_allele.items():
			norder_allele[k] = [x[0] for x in groupby(v)]

		allele_global_aln = norder_allele


		readset_ind = []
		for readname in readnameset:
			readInfo = [readname]
			for partial in allele_global_aln[readname]:
				for node in partial:
					try:
						readInfo.append((bc_variants[int(node[0])], node[1]))
					except KeyError:
						pass
			if len(readInfo) > 1:
				readset_ind.append(readInfo)
		readsets.append(readset_ind)


	return readsets	


def bc(locus_file, phase_input_files, t):
	t0 = time.clock()
	print(locus_file, phase_input_files)
	reads_dict, consec_pairs, locus_branch_mapping, bubbleCoverByRead, nodeCoveredByRead, readnames, readdict_allele, bubble_of_interest = vg_reader_multithreading(locus_file, phase_input_files, t)
	bubble_stats = open('bubble_stats.txt','w')
	total_bubble_cover = defaultdict(int)
	for ind in bubbleCoverByRead:
		for b, rs in ind.items():
			total_bubble_cover[b] += len(rs)
	for b, c in total_bubble_cover.items():
		bubble_stats.write(str(b)+'\t'+str(c)+'\n')
	bubble_stats.close()
	node_stats = open('node_stats.txt', 'w')
	for n, c in nodeCoveredByRead.items():
		node_stats.write(str(n)+'\t'+str(c)+'\n')
	node_stats.close()
	#print('reads_dict', reads_dict)
	#print('bubbleCoverByRead', bubbleCoverByRead)
	print('vg_reader took time:', time.clock() - t0)
	unitigs = find_bubble_chains(consec_pairs, nodeCoveredByRead, bubbleCoverByRead, bubble_of_interest)
	AUX = aux_contigs(unitigs, reads_dict, bubble_of_interest)
	final_ctgs = find_contigs(unitigs, AUX)
	bN = 0
	total_readsets = []
	p = Pool(t)
	processes = []
	for contig in final_ctgs:
		bubbleN = 0
		for node in contig:
			if int(node) < 0:
				bubbleN += 1
		if bubbleN >= 2:
			bN += 1
			processes.append(p.apply_async(group_readest, args = (contig, bubbleCoverByRead, readnames, readdict_allele, )))
	p.close()
	p.join()
	for proc in processes:
		readsets = proc.get()
		total_readsets.append(readsets)
	return total_readsets

