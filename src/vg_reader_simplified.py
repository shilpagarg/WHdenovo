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
	
	print('Start to read locus_file')
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
				if l.snarl.start.backward == True:
					path_in_bubble.append(tuple ((l.snarl.end.node_id,l.snarl.start.node_id)))
				else:
					path_in_bubble.append(tuple ((l.snarl.start.node_id,l.snarl.end.node_id)))
			else:
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
					trans_raw.append(l)
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
	
	locus_branch_mapping[locus_count] = per_locus
	locus_branch_mapping_raw[locus_count] = trans_raw
	het_count= 0
	for k,bubble in locus_branch_mapping.items():
		if len(bubble) > 1:
			het_count = het_count +1
	print('The number of hets:', het_count)

	reverse_mapping = defaultdict(set)
	for k, bubble in locus_branch_mapping.items():
		if bubble == []:
			continue
		for path in bubble:
			for edge in path:
				for node in edge:
					reverse_mapping[node].add(k)
	print(reverse_mapping)

	print('Reverse_mapping done')

	return reverse_mapping, locus_branch_mapping_raw, locus_branch_mapping
	# k -> bubble id, i -> path id in one bubble

def map_alignment(reverse_mapping, gam_file):

	# both simple and complex bubbles: extract reads from GAM file associated with the locus and create a sorted readset.
	# in complex bubble, set of nodes uniquely determine the path. 
	nodes_in_bubble = set()
	for i in reverse_mapping:
		nodes_in_bubble.add(i)
		#nodes_in_bubble.add(i[1])

	consec_pairs = defaultdict(set)
	reads_dict = defaultdict(list)
	rawread = defaultdict()
	bubbleCoverByRead = defaultdict(set)
	nodeCoveredByRead = defaultdict(set)
	duplicated = 0
	#TODO: consider reads with only positive score
	print('Reading alignment:', gam_file)
	with stream.open(str(gam_file), "rb") as istream:
		for data in istream:
			g = vg_pb2.Alignment()
			g.ParseFromString(data) 
			if g.score < 30:
				continue
			rawread[str(g.name) + '_' + str(g.query_position)] = g
			# hard-coded source id, mapping quality and other values.
			#sscore = g.score/len(g.sequence)

			read= [] # create read for each read alignment
			prev_tmp=[]
			prev_locus= 0
			if len(g.path.mapping) == 1:
				node = g.path.mapping[0].position.node_id
				if node in reverse_mapping:
					bubble = list(reverse_mapping[node])
					reads_dict[g.name+"_"+str(g.query_position)] = bubble
					read = bubble
					for n in bubble:
						bubbleCoverByRead[n].add(str(g.name) + '_' + str(g.query_position))
				else:
					nodeCoveredByRead[node].add(str(g.name) + '_' + str(g.query_position))
					reads_dict[g.name+"_"+str(g.query_position)].append(node)
					read.append(node)
			else:
				for i in range(0, len(g.path.mapping)):
					node = g.path.mapping[i].position.node_id
					if node in reverse_mapping:

						bubble = reverse_mapping[node]
						if len(bubble) == 1:
							try:
								if read[-1] == list(bubble)[0]:
									continue
							except IndexError:
								pass
							reads_dict[g.name+"_"+str(g.query_position)].append(list(bubble)[0])
							read.append(list(bubble)[0])
							bubbleCoverByRead[list(bubble)[0]].add(str(g.name) + '_' + str(g.query_position))
						else:
							bubble = list(bubble)
							
							if len(read) == 0:
								nextnode = g.path.mapping[i+1].position.node_id
								reads_dict[g.name+"_"+str(g.query_position)].append(bubble[0])
								reads_dict[g.name+"_"+str(g.query_position)].append(bubble[1])
								read.append(bubble[0])
								read.append(bubble[1])
								bubbleCoverByRead[bubble[0]].add(str(g.name) + '_' + str(g.query_position))
								bubbleCoverByRead[bubble[1]].add(str(g.name) + '_' + str(g.query_position))
								if list(reverse_mapping[nextnode])[0] == bubble[0]:
									reads_dict[g.name+"_"+str(g.query_position)].reverse()
									read.reverse()
							else:
								if read[-1] == bubble[0]:
									reads_dict[g.name+"_"+str(g.query_position)].append(bubble[1])
									read.append(bubble[1])
									bubbleCoverByRead[bubble[1]].add(str(g.name) + '_' + str(g.query_position))
								else:
									reads_dict[g.name+"_"+str(g.query_position)].append(bubble[0])
									read.append(bubble[0])
									bubbleCoverByRead[bubble[1]].add(str(g.name) + '_' + str(g.query_position))

					else:
						reads_dict[g.name+"_"+str(g.query_position)].append(node)
						read.append(node)
						nodeCoveredByRead[node].add(str(g.name) + '_' + str(g.query_position))
			for k in range(0,len(read)-1):
				pair1= str(read[k])+"_"+str(read[k+1]) # not taking care of reverse direction now
				pair2= str(read[k+1])+"_"+str(read[k])
				# should take of direction, not adding pairs reverse of each other
				if pair2 in consec_pairs:
					consec_pairs[pair2].add(g.name + '_' + str(g.query_position))
				else:
					consec_pairs[pair1].add(g.name + '_' + str(g.query_position))
	return reads_dict, consec_pairs, bubbleCoverByRead, rawread, nodeCoveredByRead

def vg_reader2(locus_file, gam_file):

	reverse_mapping, locus_branch_mapping_raw, locus_branch_mapping = reverse_map(locus_file)
	if type(gam_file) == str:
		print('reading alignment', gam_file)
		reads_dict, consec_pairs, bubbleCoverByRead, rawread, nodeCoveredByRead = map_alignment(reverse_mapping, gam_file)
	elif type(gam_file) == list:
		reads_dict = defaultdict()
		rawread = defaultdict()
		consec_pairs = defaultdict(set)
		bubbleCoverByRead = list()
		nodeCoveredByRead = defaultdict(set)
		for individual in gam_file:
			reads_dict_ind, consec_pairs_ind, bubbleCoverByRead_ind, rawread_ind, nodeCoveredByRead_ind = map_alignment(reverse_mapping, individual)
			reads_dict.update(reads_dict_ind)
			rawread.update(rawread_ind)
			for k, v in consec_pairs_ind.items():
				if k in consec_pairs:
					consec_pairs[k] = consec_pairs[k].union(v)
				else:
					edge = k.split('_')
					edge.reverse()
					k2 = '_'.join(edge)
					try:
						consec_pairs[k2] = consec_pairs[k2].union(v)
					except KeyError:
						consec_pairs[k2] = v
			bubbleCoverByRead.append(bubbleCoverByRead_ind)
			for k, v in nodeCoveredByRead_ind.items():
				try:
					nodeCoveredByRead[k] = nodeCoveredByRead[k].union(v)
				except KeyError:
					nodeCoveredByRead[k] = v
	return reads_dict, consec_pairs, locus_branch_mapping_raw, bubbleCoverByRead, rawread, nodeCoveredByRead
'''
############################## HARD CODED



def reverse_map2(transdata):
	
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
	for l in transdata:
		current_startsnarl = l['start']
		current_endsnarl = l['end']
		path_in_bubble = []
		if len(l['visits']) == 0:
			path_in_bubble.append((l['start'], l['end']))
		else:
			path_in_bubble.append((l['start'], l['visits'][0]))
			for i in range(len(l["visits"]) - 1):
				path_in_bubble.append((l['visits'][i], l['visits'][i+1]))
			path_in_bubble.append((l['visits'][-1], l['end']))
		if current_startsnarl == prev_startsnarl and current_endsnarl == prev_endsnarl:
			per_locus.append(path_in_bubble)
			trans_raw.append(l)
		else:
			try:
				locus_branch_mapping[locus_count] = per_locus
				locus_branch_mapping_raw[locus_count] = trans_raw
			except NameError:
				pass
			locus_count -= 1
			locus_branch_mapping[locus_count] = per_locus
			locus_branch_mapping_raw[locus_count] = trans_raw
			per_locus = []
			trans_raw = []
			per_locus.append(path_in_bubble)
			trans_raw.append(l)
		prev_startsnarl = current_startsnarl
		prev_endsnarl = current_endsnarl

	locus_branch_mapping[locus_count] = per_locus
	locus_branch_mapping_raw[locus_count] = trans_raw
	het_count= 0
	for k,bubble in locus_branch_mapping.items():
		if len(bubble) > 1:
			het_count = het_count +1
	print('The number of hets:', het_count)

	# both simple and complex bubbles: key is the values in locus_branch_mapping and value is triplet(locus, branch, alleles)
	reverse_mapping = defaultdict(set)
	for k, bubble in locus_branch_mapping.items():
		if bubble == []:
			continue
		for path in bubble:
			for edge in path:
				for node in edge:
					reverse_mapping[node].add(k)
	print(reverse_mapping)

	print('Reverse_mapping done')

	return reverse_mapping, locus_branch_mapping_raw, locus_branch_mapping
	# k -> bubble id, i -> path id in one bubble

def map_alignment2(reverse_mapping, gamdata):

	# both simple and complex bubbles: extract reads from GAM file associated with the locus and create a sorted readset.
	# in complex bubble, set of nodes uniquely determine the path. 
	nodes_in_bubble = set()
	for i in reverse_mapping:
		nodes_in_bubble.add(i)
		#nodes_in_bubble.add(i[1])

	consec_pairs = defaultdict(set)
	reads_dict = defaultdict(list)
	rawread = defaultdict()
	bubbleCoverByRead = defaultdict(set)
	nodeCoveredByRead = defaultdict(set)
	duplicated = 0
	#TODO: consider reads with only positive score
	for g in gamdata:
		rawread[str(g['name']) + '_' + str(g['query_position'])] = g
		# hard-coded source id, mapping quality and other values.
		#sscore = g.score/len(g.sequence)

		read= [] # create read for each read alignment
		prev_tmp=[]
		prev_locus= 0
		if len(g['mapping']) == 1:
			#print('DEBUG::one node raed', g.name)
			node = g['mapping'][0]
			if node in reverse_mapping:
				bubble = list(reverse_mapping[node])
				reads_dict[g['name']+"_"+str(g['query_position'])] = bubble
				read = bubble
				for n in bubble:
					bubbleCoverByRead[n].add(str(g['name']) + '_' + str(g['query_position']))
			else:
				bubbleCoverByRead[node].add(str(g['name']) + '_' + str(g['query_position']))
				reads_dict[g['name']+"_"+str(g['query_position'])].append(node)
				read.append(node)
			#print('DEBUG::one node read', read)
		else:
			#print('DEBUG::multi-node read', g.name)
			for i in range(0, len(g['mapping'])):
				node = g['mapping'][i]
				if node in reverse_mapping:

					bubble = reverse_mapping[node]

					if len(bubble) == 1:
						try:
							if read[-1] == list(bubble)[0]:
								continue
						except IndexError:
							pass
						reads_dict[str(g['name']) + '_' + str(g['query_position'])].append(list(bubble)[0])
						read.append(list(bubble)[0])
						nodeCoveredByRead[list(bubble)[0]].add(str(g['name']) + '_' + str(g['query_position']))
					else:
						bubble = list(bubble)
						
						if len(read) == 0:
							nextnode = g['mapping'][i+1]
							reads_dict[str(g['name']) + '_' + str(g['query_position'])].append(bubble[0])
							reads_dict[str(g['name']) + '_' + str(g['query_position'])].append(bubble[1])
							read.append(bubble[0])
							read.append(bubble[1])
							bubbleCoverByRead[bubble[0]].add(str(g['name']) + '_' + str(g['query_position']))
							bubbleCoverByRead[bubble[1]].add(str(g['name']) + '_' + str(g['query_position']))
							if list(reverse_mapping[nextnode])[0] == bubble[0]:
								reads_dict[str(g['name']) + '_' + str(g['query_position'])].reverse()
								read.reverse()
						else:
							if read[-1] == bubble[0]:
								reads_dict[str(g['name']) + '_' + str(g['query_position'])].append(bubble[1])
								read.append(bubble[1])
								bubbleCoverByRead[bubble[1]].add(str(g['name']) + '_' + str(g['query_position']))
							else:
								reads_dict[str(g['name']) + '_' + str(g['query_position'])].append(bubble[0])
								read.append(bubble[0])
								bubbleCoverByRead[bubble[1]].add(str(g['name']) + '_' + str(g['query_position']))

				else:
					reads_dict[str(g['name']) + '_' + str(g['query_position'])].append(node)
					read.append(node)
					bubbleCoverByRead[node].add(str(g['name']) + '_' + str(g['query_position']))
					nodeCoveredByRead[node].add(str(g['name']) + '_' + str(g['query_position']))

		if len(read) == 0:
			print('DEBUG::len(read) == 0:', g.name)
		for k in range(0,len(read)-1):
			pair1= str(read[k])+"_"+str(read[k+1]) # not taking care of reverse direction now
			pair2= str(read[k+1])+"_"+str(read[k])
			# should take of direction, not adding pairs reverse of each other
			if pair2 in consec_pairs:
				consec_pairs[pair2].add(str(g['name']) + '_' + str(g['query_position']))
			else:
				consec_pairs[pair1].add(str(g['name']) + '_' + str(g['query_position']))

	return reads_dict, consec_pairs, bubbleCoverByRead, rawread, nodeCoveredByRead


def vg_reader_hard():
	start = 'start'
	end = 'end'
	transdata = []
	transdata.append({'visits':[15], 'start':14, 'end': 17})
	transdata.append({'visits':[16], 'start':14, 'end': 17})

	name = 'name'
	query_position = 'query_position'
	mapping = 'mapping'
	gamdata = []
	gamdata.append({name: 'r1', query_position: 0, mapping: [1, 2, 3, 4]})
	gamdata.append({name: 'r2', query_position: 0, mapping: [4, 3, 5, 6, 7]})
	gamdata.append({name: 'r3', query_position: 0, mapping: [10, 9, 8, 7, 6]})
	gamdata.append({name: 'r2', query_position: 0, mapping: [8, 12, 13, 14, 16, 17, 18]})
	gamdata.append({name: 'r4', query_position: 0, mapping: [5, 6]})
	gamdata.append({name: 'r4', query_position: 3, mapping: [8, 9]})
	gamdata.append({name: 'r5', query_position: 1, mapping: [14, 15, 17, 18, 19,20, 21, 22]})

	reverse_mapping, locus_branch_mapping_raw, locus_branch_mapping = reverse_map2(transdata)
	print('DEBUG::locus_branch_mapping', locus_branch_mapping)
	
	reads_dict, consec_pairs, bubbleCoverByRead, rawread, nodeCoveredByRead = map_alignment2(reverse_mapping, gamdata)
	print(reads_dict)
	print(consec_pairs)
	return reads_dict, consec_pairs, locus_branch_mapping_raw, bubbleCoverByRead, rawread, nodeCoveredByRead
############################## HARD CODED
'''
