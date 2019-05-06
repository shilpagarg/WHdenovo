from collections import defaultdict, OrderedDict
import stream
import vg_pb2

def mergePath(tempPath, path_in_bubble, insideBack, pathBack, local_path_back):
	out = []
	if not insideBack:
		insideBack = -1
	else:
		insideBack = 1
	if not pathBack:
		pathBack = -1
	else:
		pathBack = 1
	if insideBack * pathBack * local_path_back == -1:
		newpath = []
		path_in_bubble.reverse()
		for n in path_in_bubble:
			newpath.append((n[1], n[0]))
		path_in_bubble = newpath.copy()

	for i in range(len(tempPath)):
		if 0 not in tempPath[i]:
			out.append(tempPath[i])
		else:
			break
	if pathBack == -1:
		out.append((tempPath[i][0], path_in_bubble[0][0]))
	for j in path_in_bubble:
		out.append(j)
	if pathBack == 1:
		out.append((j[1], tempPath[i+1][1]))
	for k in range(i + 2, len(tempPath)):
		out.append(tempPath[k])
	return out

def reverse_map(locus_file):

	locus_count = 0
	per_locus = []
	trans_raw = []
	prev_startsnarl = 0
	prev_endsnarl = 0
	locus_branch_mapping = OrderedDict()
	locus_branch_mapping_raw = OrderedDict()
	prev_startsnarl_orientation = -1
	prev_endsnarl_orientation = -1
	insidebubble = 0
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
			path_in_bubble = []
			hasInBubble = False


			if len(l.visits) == 0:
				#TODO: for now, assumed, all nodes in path are either forward or backward
				if l.snarl.start.backward == True:
					path_in_bubble.append(tuple ((l.snarl.end.node_id,l.snarl.start.node_id)))
				else:
					path_in_bubble.append(tuple ((l.snarl.start.node_id,l.snarl.end.node_id)))
			else:
				#TODO: for now, assumed, all nodes in path are either forward or backward
				if (l.snarl.start.backward == True and l.snarl.end.backward != True) or (l.snarl.start.backward != True and l.snarl.end.backward == True):
					path_in_bubble.append(tuple ((l.snarl.end.node_id, l.visits[-1].node_id)))
					local_path_back = -1
					for i in range(len(l.visits)):
						if l.visits[i].snarl.start.node_id != 0:
							pathBack = True
							if l.visits[i].backward:
								insideBack = True
							else:
								insideBack = False
							insidebubble = 1
							hasInBubble = True
						if i == len(l.visits) - 1:
							break
						path_in_bubble.append(tuple((l.visits[-1 - i].node_id, l.visits[-2 - i].node_id)))
					path_in_bubble.append(tuple ((l.visits[0].node_id,l.snarl.start.node_id)))
				else:
					local_path_back = 1
					path_in_bubble.append(tuple ((l.snarl.start.node_id,l.visits[0].node_id)))
					for i in range(len(l.visits)):
						if l.visits[i].snarl.start.node_id != 0:
							pathBack = False
							if l.visits[i].backward:
								insideBack = True
							else:
								insideBack = False
							insidebubble = 1
							hasInBubble = True
						if i == len(l.visits) - 1:
							break
						path_in_bubble.append(tuple((l.visits[i].node_id, l.visits[i + 1].node_id)))
						
					path_in_bubble.append(tuple ((l.visits[-1].node_id, l.snarl.end.node_id))) 
						
			if hasInBubble:
				tempPath = path_in_bubble.copy()

				if current_startsnarl == prev_startsnarl and current_endsnarl == prev_endsnarl and current_endsnarl_orientation == prev_endsnarl_orientation and prev_startsnarl_orientation == current_startsnarl_orientation:
					pass
				else:
					try:
						locus_branch_mapping[locus_count] = per_locus
						locus_branch_mapping_raw[locus_count] = trans_raw
					except NameError:
						pass
					locus_count -= 1
					per_locus = []
					trans_raw = []
			else:
				if current_startsnarl == prev_startsnarl and current_endsnarl == prev_endsnarl and current_endsnarl_orientation == prev_endsnarl_orientation and prev_startsnarl_orientation == current_startsnarl_orientation:
					if insidebubble == 2:
						path_in_bubble = mergePath(tempPath, path_in_bubble, insideBack, pathBack, local_path_back)
						per_locus.append(path_in_bubble)
						trans_raw.append(l)
						insidebubble = 0
						insideBack = False
						pathBack = False
					else:
						per_locus.append(path_in_bubble)
						trans_raw.append(l)
				else:

					if insidebubble == 1:
						insidebubble = 2
						path_in_bubble = mergePath(tempPath, path_in_bubble, insideBack, pathBack, local_path_back)
						per_locus.append(path_in_bubble)
						trans_raw.append(l)
					else:
						try:
							locus_branch_mapping[locus_count] = per_locus
							locus_branch_mapping_raw[locus_count] = trans_raw
						except NameError:
							pass
						locus_count -= 1
						per_locus = []
						trans_raw = []
						per_locus.append(path_in_bubble)
						trans_raw.append(l)

			prev_startsnarl = current_startsnarl
			prev_startsnarl_orientation = current_startsnarl_orientation
			prev_endsnarl = current_endsnarl
			prev_endsnarl_orientation = current_endsnarl_orientation
	
	print('The number of hets:')
	het_count= 0
	for k,bubble in locus_branch_mapping.items():
		if len(bubble) > 1:
			het_count = het_count +1
	print(het_count)
	# keep branch of paths in each bubble.
	alleles_per_pos = defaultdict()
	for k, bubble in locus_branch_mapping.items():
		alleles_per_pos[k] = len(bubble)

	# both simple and complex bubbles: key is the values in locus_branch_mapping and value is triplet(locus, branch, alleles)
	reverse_mapping = defaultdict(list)
	for k, bubble in locus_branch_mapping.items():
		if len(bubble) > 1: # more than one branch
			for i, path in enumerate(bubble):
				if len(path) > 0:
					for edge in path:
						reverse_mapping[edge].append([k, i, len(path), len(bubble)]) # in complex bubbles, a node can map to multiple branches.


	return reverse_mapping, locus_branch_mapping_raw, locus_branch_mapping
	# k -> bubble id, i -> path id in one bubble
	print('reverse_mapping')
	print(reverse_mapping)

def map_alignment(reverse_mapping, gam_file):

	# both simple and complex bubbles: extract reads from GAM file associated with the locus and create a sorted readset.
	# in complex bubble, set of nodes uniquely determine the path. 
	nodes_in_bubble = set()
	for i in reverse_mapping:
		nodes_in_bubble.add(i[0])
		nodes_in_bubble.add(i[1])

	consec_pairs = defaultdict(set)
	triples = defaultdict(set)
	reads_dict = defaultdict(list)
	rawread = defaultdict()
	bubbleCoverByRead = defaultdict(set)
	duplicated = 0
	#TODO: consider reads with only positive score.
	with stream.open(str(gam_file), "rb") as istream:
		for data in istream:
			g = vg_pb2.Alignment()
			g.ParseFromString(data) 
			rawread[str(g.name) + '_' + str(g.query_position)] = g
			# hard-coded source id, mapping quality and other values.
			val1 = True
			val2 = False

			count1 =0
			count2=0
			score = g.score/len(g.sequence)

			read= [] # create read for each read alignment
			prev_tmp=[]
			prev_locus= 0
			for i in range(0,len(g.path.mapping)-1):
				edge1 = tuple((int(g.path.mapping[i].position.node_id), int(g.path.mapping[i+1].position.node_id))) # go over nodes in a mapping
				edge2 = tuple((int(g.path.mapping[i+1].position.node_id), int(g.path.mapping[i].position.node_id))) # go over nodes in a mapping
				if edge1 in reverse_mapping or edge2 in reverse_mapping: # handle start and sink node.
					if edge1 in reverse_mapping:
						#qualities = [10]* reverse_mapping[edge1][0][2] # reverse_mapping[edge1][0][2] -> the number of paths of the bubble this edge belongs to.
						node_inf = [tuple(j[0:3]) for j in reverse_mapping[edge1]] # consider (locus, branch)
					else:
						node_inf = [tuple(j[0:3]) for j in reverse_mapping[edge2]]
						# node_inf = [(bubble_id, path_id_in_bubble, len_of_this_path), ( , , ), ...]
					tmp = node_inf.copy()
					if prev_locus != tmp[0][0]:
						# entering a new bubble
						prev_tmp = tmp.copy()
						prev_locus = tmp[0][0]
						len_in_path = 1
					else:
						len_in_path += 1
					interset_tmp = list(set(tmp).intersection(set(prev_tmp)))
					if len(interset_tmp) == 1 and interset_tmp[0][2] == len_in_path:
						reads_dict[g.name+"_"+str(g.query_position)].append(interset_tmp[0][0])
						read.append(interset_tmp[0][0])
						
						bubbleCoverByRead[interset_tmp[0][0]].add(str(g.name) + '_' + str(g.query_position))
					prev_tmp = interset_tmp.copy()
				else:
					if i == 0:
						if int(g.path.mapping[i].position.node_id) not in nodes_in_bubble:
							read.append(int(g.path.mapping[i].position.node_id))
							reads_dict[g.name+"_"+str(g.query_position)].append(int(g.path.mapping[i].position.node_id))

					if int(g.path.mapping[i + 1].position.node_id) not in nodes_in_bubble:
						read.append(int(g.path.mapping[i + 1].position.node_id))
						reads_dict[g.name+"_"+str(g.query_position)].append(int(g.path.mapping[i + 1].position.node_id))
			# for every pair of bubbles or bubble-node			
			for k in range(0,len(read)-1):
				pair1= str(read[k])+"_"+str(read[k+1]) # not taking care of reverse direction now
				pair2= str(read[k+1])+"_"+str(read[k])
				# should take of direction, not adding pairs reverse of each other
				if pair2 in consec_pairs:
					consec_pairs[pair2].add(g.name + '_' + str(g.query_position))
				else:
					consec_pairs[pair1].add(g.name + '_' + str(g.query_position))
			for k in range(len(read) - 2):
				triple1 = (read[k], read[k+1], read[k+2])
				triple2 = (read[k+2], read[k+1], read[k])
				if triple1 in triples:
					triples[triple1].add(g.name)
				else:
					triples[triple2].add(g.name)


	return reads_dict, consec_pairs, triples, bubbleCoverByRead, rawread

def vg_reader2(locus_file, gam_file):

	reverse_mapping, locus_branch_mapping_raw, locus_branch_mapping = reverse_map(locus_file)
	if type(gam_file) == str:
		reads_dict, consec_pairs, triples, bubbleCoverByRead, rawread = map_alignment(reverse_mapping, gam_file)
	elif type(gam_file) == list:
		reads_dict = defaultdict()
		rawread = defaultdict()
		consec_pairs = defaultdict(set)
		triples = defaultdict(set)
		bubbleCoverByRead = list()
		for individual in gam_file:
			reads_dict_ind, consec_pairs_ind, triples_ind, bubbleCoverByRead_ind, rawread_ind = map_alignment(reverse_mapping, individual)
			reads_dict.update(reads_dict_ind)
			rawread.update(rawread_ind)
			for k, v in consec_pairs_ind.items():
				try:
					consec_pairs[k] = consec_pairs[k].union(v)
				except KeyError:
					consec_pairs[k] = v
			for k, v in triples_ind.items():
				try:
					triples[k] = triples[k].union(v)
				except KeyError:
					triples[k] = v
			bubbleCoverByRead.append(bubbleCoverByRead_ind)
	return reads_dict, consec_pairs, locus_branch_mapping_raw, bubbleCoverByRead, rawread

