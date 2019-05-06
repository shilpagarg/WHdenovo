import sys
import stream
import vg_pb2
from collections import Counter
from collections import defaultdict, OrderedDict
import networkx as nx
from itertools import groupby
from vg_reader_4v import vg_reader2

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
	consec_pairs_final=defaultdict(set)

	# pairs to trust
	for i,j in consec_pairs.items():
		#if len(j) > 0 and len(j) < 45:
		# if len(j) > 8 and len(j) < 200:
		if len(j) > 5 and len(j) < 25:
			consec_pairs_final[i] = j

	print('the total number of edges %d after removing errorenous ones.' %len(consec_pairs_final))
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

	count = 0
	for i,j in consec_pairs_tmp.items():
		if len(j) > 2:
			count+=1
			#print(i,j)
			for k in list(j)[1:]:  # instead of for k in j, randomly maintain one edge here.
				if str(i)+"_"+str(k) in consec_pairs_final:
					del consec_pairs_final[str(i)+"_"+str(k)]
				if str(k)+"_"+str(i) in consec_pairs_final:
					del consec_pairs_final[str(k)+"_"+str(i)]

	print('the total number of branching points %d.' %count)
	print('the total number of good pairs or edges %d.' %len(consec_pairs_final))

	# Return all distinct simple paths ending at "targetNode", continuing
	# from "currentPath". "usedNodes" is useful so we can quickly skip
	# nodes we have already added to "currentPath". When a new solution path
	# is found, append it to "answerPaths" and return it.
	#	


	nodeToNodess= defaultdict(set)
	nodeToNodes = defaultdict(list)
	for i,j in consec_pairs_final.items():
		if i.split("_")[0] != i.split("_")[1]:
			nodeToNodess[i.split("_")[0]].add(i.split("_")[1])
			nodeToNodess[i.split("_")[1]].add(i.split("_")[0])

	for i,j in nodeToNodess.items():
		nodeToNodes[i] = list(j)

	for i,j in nodeToNodes.items():
		if len(nodeToNodes[i]) ==1:
			start_or_end.add(i)
			
	print('the total number of start or end nodes %d.' %len(start_or_end))
	nx_nodeToNodes = nx.MultiDiGraph(nodeToNodes)

	unitigs = []
	hashse = {}
	for i in start_or_end:
		hashse[i]=0

	for i in start_or_end:
		if hashse[i]==0:
			unitigs.append([])
			for j in dfs(nodeToNodess, str(i)):
				unitigs[-1].append(j)
				if j in hashse:
					hashse[j]=1
	return unitigs

def aux_contigs(unitigs, reads_dict):
	unitigs_ends = {}
	unitigs_starts = {}
	bubble_count = 0
	bubble_dict = {}

	for var in unitigs:
		if len(var) < 2 & len(var) != 0:
			unitigs_ends[var[0]] = int(var[0])
			unitigs_starts[var[0]] = int(var[0])
		else:
		    unitigs_ends[var[0]] = int(var[-1])
		    unitigs_starts[var[-1]] = int(var[0])
		for i in range(len(var)):
			if int(var[i]) < 300:
				bubble_count += 1
				if var[i] in bubble_dict.keys():
					bubble_dict[var[i]] += 1
				else:
					bubble_dict[var[i]] = 1


	print("UNITIG ENDS \n")

	unitigs_enddict = defaultdict(int)

	for i in unitigs_ends.keys():
		unitigs_enddict[i] = 0
	for i in unitigs_starts.keys():
		unitigs_enddict[i] = 0

	nodes_list = set()
	dummy_list = ['0']*10
	orderalignment = defaultdict(list)
	orderalignment = defaultdict(lambda: [-1]*10, orderalignment)

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

	aux_unitigs = []
	bubble_dict2 = dict.fromkeys(bubble_dict.keys(), 0)
	for k,v in new_orderalignment.items():
		new_g = []
		count =0
		for i in range(0,len(v)):
			for j in range(0,len(v[i])):
				if str(v[i][j]) in bubble_dict2.keys():
					if bubble_dict2[str(v[i][j])] == 0:
						bubble_dict2[str(v[i][j])] = 1
				if str(v[i][j]) in unitigs_ends.keys():				
					if unitigs_enddict[v[i][j]] == 0:
						new_g.append(unitigs_starts[str(unitigs_ends[str(v[i][j])])])
						new_g.append(unitigs_ends[str(v[i][j])])
						unitigs_enddict[v[i][j]] = 1
						unitigs_enddict[unitigs_ends[str(v[i][j])]] = 1
				if str(v[i][j]) in unitigs_starts.keys():				
					if unitigs_enddict[v[i][j]] == 0:
						new_g.append(unitigs_ends[str(unitigs_starts[str(v[i][j])])])
						new_g.append(unitigs_starts[str(v[i][j])])
						unitigs_enddict[v[i][j]] = 1
						unitigs_enddict[unitigs_starts[str(v[i][j])]] = 1
		if len(new_g) > 0:
			aux_unitigs.append(new_g)

	return aux_unitigs 

def find_contigs(unitigs, aux_unitigs):
	consec_pairs_final = defaultdict(set)
	consec_pairs = defaultdict(set)

	count=0   
	nodeToNodess= defaultdict(set)
	nodeToNodes = defaultdict(list)

	for var in aux_unitigs:
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
		if int(key) < 300:
			count += 1
			bubble1.append(key)

	start_or_end = set()
	for i,j in nodeToNodes.items():
		if len(nodeToNodes[i]) ==1:
			start_or_end.add(i)
	print('the total number of start or end nodes %d.' %len(start_or_end))
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

	return final_ctgs

def output_trans_gam(locus_branch_mapping_raw, final_ctgs, prefix, bubbleCoverByRead, rawread):
	n_contig = 0
	for contig in final_ctgs:
		n_bubble = 0
		for i in contig:
			if int(i) < 0:
				n_bubble += 1
		if n_bubble < 2: # Filter on "at least having 2 bubbles in the output trans file"
			continue
		if len(contig) > 2: # HERE WE CAN ALSO LOST BUBBLES
			n_contig += 1
			ostream = stream.open(prefix + '_contigs_%d.trans' % n_contig, 'wb')
			if type(bubbleCoverByRead) == defaultdict:
				gamstream = stream.open(prefix + '_contigs_%d.gam' % n_contig, 'wb')
				readset = set()
			if type(bubbleCoverByRead) == list:
				gam0 = stream.open(prefix + '_contigs_%d_0.gam' % n_contig, 'wb')
				gam1 = stream.open(prefix + '_contigs_%d_1.gam' % n_contig, 'wb')
				gam2 = stream.open(prefix + '_contigs_%d_2.gam' % n_contig, 'wb')
				readset0 = set()
				readset1 = set()
				readset2 = set()
			for i in range(len(contig)):
				if type(contig[i]) == tuple:
					node = int(contig[i][0])
				else:
					node = int(contig[i])
				if node >= 0:
					continue
				try:
					for k in locus_branch_mapping_raw[node]:
						ostream.write(k)
				except KeyError:
					pass
				if type(bubbleCoverByRead) == defaultdict:
					try:
						readset = readset.union(bubbleCoverByRead[node])
					except KeyError:
						pass
				if type(bubbleCoverByRead) == list:
					try:
						readset0 = readset0.union(bubbleCoverByRead[0][node])
					except KeyError:
						pass
					try:
						readset1 = readset1.union(bubbleCoverByRead[1][node])
					except KeyError:
						pass
					try:
						readset2 = readset2.union(bubbleCoverByRead[2][node])
					except KeyError:
						pass
			ostream.close()
			try:
				for i in readset:
					gamstream.write(rawread[i])
				gamstream.close()
			except NameError:
				for i in readset0:
					gam0.write(rawread[i])
				gam0.close()
				for i in readset1:
					gam1.write(rawread[i])
				gam1.close()
				for i in readset2:
					gam2.write(rawread[i])
				gam2.close()

def main():
	locus_file = sys.argv[1]
	if len(sys.argv) == 6:
		gam_file = [sys.argv[2], sys.argv[3], sys.argv[4]]
		prefix = sys.argv[5]
	elif len(sys.argv) == 4:
		gam_file = sys.argv[2]
		prefix = sys.argv[3]
	else:
		print('Input incorrectly. Should be either (for individual):'
			  'python bubble_chain_4scripts_version.py asm1.trans aln.gam <output prefix>'
			  'Or (for trio):'
			  'python bubble_chain_4scripts_version.py asm1.trans aln0.gam aln1.gam aln2.gam <output prefix>')
	reads_dict, consec_pairs, locus_branch_mapping_raw, bubbleCoverByRead, rawread = vg_reader2(locus_file, gam_file)
	unitigs = find_bubble_chains(consec_pairs)
	aux_unitigs = aux_contigs(unitigs, reads_dict)
	final_ctgs = find_contigs(unitigs, aux_unitigs)
	output_trans_gam(locus_branch_mapping_raw, final_ctgs, prefix, bubbleCoverByRead, rawread)

main()
