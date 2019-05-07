import sys
import stream
import logging
import vg_pb2 
from collections import Counter
from collections import defaultdict
import collections
from collections import OrderedDict
from itertools import groupby

file_input = sys.argv[1]
locus_input = sys.argv[4]
unitigs = sys.argv[2]


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

	print('The number of hets:')
	het_count= 0
	for k,v in locus_branch_mapping.items():
		if len(v) >1:
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
			#       continue
			read= [] # create read for each read alignment

			prev_tmp=[]
			prev_locus= -1
			locus = -1

			for i in range(0,len(g.path.mapping)-1):
			#for i in g.path.mapping: # go over the mapping in a read
				# TODO: check for forward or reverse strand, we may not need it for DAG.
				edge1 = tuple((int(g.path.mapping[i].position.name), int(g.path.mapping[i+1].position.name))) # go over nodes in a mapping
				edge2 = tuple((int(g.path.mapping[i+1].position.name), int(g.path.mapping[i].position.name))) # go over nodes in a mapping
				if edge1 in reverse_mapping or edge2 in reverse_mapping: # handle start and sink node.
					if edge1 in reverse_mapping:
						qualities = [10]* reverse_mapping[edge1][0][2]
						node_inf= [tuple(i[0:2]) for i in reverse_mapping[edge1]] # consider (locus, branch)
					else:
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
						else:
							next_edge1 = tuple((int(g.path.mapping[i+1].position.name), int(g.path.mapping[i+2].position.name)))
							next_edge2 = tuple((int(g.path.mapping[i+2].position.name), int(g.path.mapping[i+1].position.name)))

							if next_edge1 not in reverse_mapping and next_edge2 not in reverse_mapping:
								#read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)    
								reads_dict[g.name+"_"+str(g.query_position)].append(interset_tmp[0][0])
								read.append(interset_tmp[0][0])
						locus= interset_tmp[0][0]
				else:
					read.append(int(g.path.mapping[i].position.name))
					read.append(int(g.path.mapping[i+1].position.name))
					reads_dict[g.name+"_"+str(g.query_position)].append(int(g.path.mapping[i].position.name))
					reads_dict[g.name+"_"+str(g.query_position)].append(int(g.path.mapping[i+1].position.name))

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


reads_dict, consec_pairs = vg_reader(locus_input, file_input)
print("READS DICT", reads_dict)
def reverse_complement(seq):
	seq_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'T', 'c': 'G', 'g': 'C', 't': 'A'}
	return "".join([seq_dict[base] for base in reversed(seq)])

unitigs_ends = {}
unitigs_starts = {}
#print(unitigs)
bubble_count = 0
bubble_dict = {}
with open(unitigs) as f:
        for line in f:
                var = line.rstrip().split(",")
                #print (var)
                # take care of last comma
                if len(var) < 3 & len(var) != 0:
                	unitigs_ends[var[0]] = int(var[0])
                	unitigs_starts[var[0]] = int(var[0])
                else:
                    unitigs_ends[var[0]] = int(var[-2])
                    unitigs_starts[var[-2]] = int(var[0])
                for i in range(len(var)-1):
                	if int(var[i]) < 300:
                		bubble_count += 1
                		if var[i] in bubble_dict.keys():
                			bubble_dict[var[i]] += 1
                		else:
                			bubble_dict[var[i]] = 1
print("UNITIG ENDS \n")
print(len(unitigs_ends))
print(len(unitigs_starts))
print(bubble_count)
print(bubble_dict)
unitigs_enddict = defaultdict(int)

for i in unitigs_ends.keys():
	unitigs_enddict[i] = 0
for i in unitigs_starts.keys():
	unitigs_enddict[i] = 0

#print(unitigs_enddict)

nodes_list = set()
dummy_list = ['0']*10
orderalignment = defaultdict(list)
orderalignment = defaultdict(lambda: [-1]*10, orderalignment)

for i,j in reads_dict.items():
	canu_name = '_'.join(i.split("_")[:-1])
	canu_chunk_num = int(i.split("_")[-1])
	orderalignment[canu_name].insert(canu_chunk_num, j)

#print("ORDERALIGNMENT", orderalignment)

new_orderalignment = defaultdict(list)
for k,v in orderalignment.items():
	new_orderalignment[k] = [x for x in v if x != -1]

norder = defaultdict(list)
for k,v in new_orderalignment.items():
	#print("VVVVVVVVV", k, v)
	norder[k] = [x[0] for x in groupby(v)]

new_orderalignment = norder

#print("NORDER \n", norder)
#print("LENGTH", len(new_orderalignment))

out = open(sys.argv[3], 'w')
bubble_dict2 = dict.fromkeys(bubble_dict.keys(), 0)
for k,v in new_orderalignment.items():
	new_g = []
	count =0
	for i in range(0,len(v)):
		#print(v[i])
		#new_g.name = k
		#tmp = v[i].name.split("_")[0]
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
	for p in new_g:
		out.write(str(p) + ",")
	if len(new_g) > 0:
		out.write("\n")
print(sum(bubble_dict2.values()), '\n', bubble_dict2)
