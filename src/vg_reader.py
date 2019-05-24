import sys
import subprocess
from multiprocessing import Pool, Queue
import time
import math
import json
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
	#trans_raw = []
	prev_startsnarl = 0
	prev_endsnarl = 0
	locus_branch_mapping = OrderedDict()
	#locus_branch_mapping_raw = OrderedDict()
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
					#trans_raw.append(l)
					pass
				else:
					try:
						locus_branch_mapping[locus_count] = per_locus
						#locus_branch_mapping_raw[locus_count] = trans_raw
					except NameError:
						pass
					locus_count -= 1
					per_locus = []
					#trans_raw = []
			else:
				if current_startsnarl == prev_startsnarl and current_endsnarl == prev_endsnarl and current_endsnarl_orientation == prev_endsnarl_orientation and prev_startsnarl_orientation == current_startsnarl_orientation:
					if insidebubble == 2:
						path_in_bubble = mergePath(tempPath, path_in_bubble, insideBack, pathBack, local_path_back)
						per_locus.append(path_in_bubble)
						#trans_raw.append(l)
						insidebubble = 0
						insideBack = False
						pathBack = False
					else:
						per_locus.append(path_in_bubble)
						#trans_raw.append(l)
				else:
					if insidebubble == 1:
						insidebubble = 2
						path_in_bubble = mergePath(tempPath, path_in_bubble, insideBack, pathBack, local_path_back)
						per_locus.append(path_in_bubble)
						#trans_raw.append(l)
					else:
						try:
							locus_branch_mapping[locus_count] = per_locus
							#locus_branch_mapping_raw[locus_count] = trans_raw
						except NameError:
							pass
						locus_count -= 1
						per_locus = []
						per_locus.append(path_in_bubble)

			prev_startsnarl = current_startsnarl
			prev_startsnarl_orientation = current_startsnarl_orientation
			prev_endsnarl = current_endsnarl
			prev_endsnarl_orientation = current_endsnarl_orientation
	
	locus_branch_mapping[locus_count] = per_locus
	#locus_branch_mapping_raw[locus_count] = trans_raw
	het_count= 0
	for k,bubble in locus_branch_mapping.items():
		if len(bubble) > 1:
			het_count = het_count +1
	print('The number of hets:', het_count)

	reverse_mapping = defaultdict(set)
	allele_reverse_mapping = defaultdict(list)
	for k, bubble in locus_branch_mapping.items():
		if bubble == []:
			continue
		for path in bubble:
			for edge in path:
				for node in edge:
					reverse_mapping[node].add(k)
		for i, path in enumerate(bubble):
			if len(path) > 0:
				for edge in path:
					allele_reverse_mapping[edge].append([k, i, len(path), len(bubble)])
	#print('DEBUG::locus_branch_mapping',locus_branch_mapping)
	print('Reverse_mapping done')

	return reverse_mapping, locus_branch_mapping, allele_reverse_mapping

def alleleDetect(contents, allele_reverse_mapping, chunk_id, chunkSize):

	readdict_allele = defaultdict(list)
	for p in range(len(contents)):
		partial = contents[p]
		partial_line = chunk_id * chunkSize + p + 1 # TODO: This is a 1-based row number in original json. See if better to be 0-based.
		g = json.loads(partial.rstrip())
		try:
			name = g['name']
		except KeyError:
			return None	
		try:
			qp = g['query_position']
		except KeyError:
			qp = 0
		try:
			score = g['score']
		except KeyError:
			score = 0
		mapping = g['path']['mapping']
		rn_qp = name + '_' + str(qp)

		#score = g.score/len(g.sequence)
		#if score > 0.2:
		#   continue
		prev_tmp = []
		prev_locus = -1
		n_variant = 0
		added = False

		for i in range(len(mapping) - 1):
		#for i in g.path.mapping: # go over the mapping in a read
		# TODO: check for forward or reverse strand, we may not need it for DAG.

			edge1 = (mapping[i]['position']['node_id'], mapping[i+1]['position']['node_id']) # go over nodes in a mapping
			edge2 = (mapping[i+1]['position']['node_id'], mapping[i]['position']['node_id']) # go over nodes in a mapping

			if edge1 in allele_reverse_mapping or edge2 in allele_reverse_mapping: # handle start and sink node.
				if edge1 in allele_reverse_mapping:   
					#qualities = [10]* reverse_mapping[edge1][0][2]
					qualitie = 1 
					node_inf = [tuple(i[0:3]) for i in allele_reverse_mapping[edge1]] # consider (locus, branch)
				else:
					# qualities = [10]* reverse_mapping[edge2][0][2]
					qualities = 1 
					node_inf = [tuple(i[0:3]) for i in allele_reverse_mapping[edge2]]
				tmp = node_inf.copy()
				if prev_locus != tmp[0][0]:
					added = False
					prev_tmp = tmp.copy()
					prev_locus = tmp[0][0]
					#len_in_path = 1
				#else:
					#len_in_path += 1
				if added:
					continue
				interset_tmp = list(set(tmp).intersection(set(prev_tmp)))
				if len(interset_tmp) == 1:# and interset_tmp[0][2] == len_in_path: # for complicated bubbles, but with Top-k paths. combination of some nodes uniquely determine branch.
					qualities = 1 
					readdict_allele[rn_qp].append((interset_tmp[0][0], interset_tmp[0][1]))
					#read.add_variant(interset_tmp[0][0], interset_tmp[0][1], qualities)
					added = True
	return readdict_allele


def AlignmentParse(contents, reverse_mapping, chunk_id, chunkSize):

	consec_pairs = set()
	#consec_pairs = defaultdict(set)
	reads_dict = defaultdict(list)
	readnames = []
	bubbleCoverByRead = defaultdict(set)
	nodeCoveredByRead = defaultdict(int)
	for p in range(len(contents)):
		partial = contents[p]
		partial_line = chunk_id * chunkSize + p + 1 # TODO: This is a 1-based row number in original json. See if better to be 0-based.
		g = json.loads(partial.rstrip())
		try:
			name = g['name']
		except KeyError:
			return None
		try:
			qp = g['query_position']
		except KeyError:
			qp = 0
		try:
			score = g['score']
		except KeyError:
			score = 0
		readnames.append(name)
		mapping = g['path']['mapping']
		rn_qp = name + '_' + str(qp)
		hasDup = False
		if len(mapping) == 1:
			node = mapping[0]['position']['node_id']
			if node in reverse_mapping:
				bubble = list(reverse_mapping[node])
				reads_dict[rn_qp] = bubble
				for n in bubble:
					bubbleCoverByRead[n].add(partial_line)
			else:
				nodeCoveredByRead[node] += 1
				#nodeCoveredByRead[node].add(rn_qp)
				reads_dict[rn_qp].append(node)
		else:
			seenNode = set()
			for i in range(len(mapping)):
				node = mapping[i]['position']['node_id']
				if node in reverse_mapping:
					bubble = reverse_mapping[node]
					if len(bubble) == 1:
						try:
							if reads_dict[rn_qp][-1] == list(bubble)[0]:
								continue
						except IndexError:
							pass
						reads_dict[rn_qp].append(list(bubble)[0])
						if list(bubble)[0] in seenNode:
							hasDup = True
							#break
						seenNode.add(list(bubble)[0])
						bubbleCoverByRead[list(bubble)[0]].add(partial_line)
					else:
						bubble = list(bubble)
						if len(reads_dict[rn_qp]) == 0:
							nextnode = mapping[i+1]['position']['node_id']
							if bubble[0] in seenNode or bubble[1] in seenNode:
								hasDup = True
								#break
							reads_dict[rn_qp].append(bubble[0])
							seenNode.add(bubble[0])
							reads_dict[rn_qp].append(bubble[1])
							seenNode.add(bubble[1])
							bubbleCoverByRead[bubble[0]].add(partial_line)
							bubbleCoverByRead[bubble[1]].add(partial_line)
							if list(reverse_mapping[nextnode])[0] == bubble[0]:
								reads_dict[rn_qp].reverse()
						else:
							if reads_dict[rn_qp][-1] == bubble[0]:
								if bubble[1] in seenNode:
									hasDup = True
									#break
								reads_dict[rn_qp].append(bubble[1])
								seenNode.add(bubble[1])
								bubbleCoverByRead[bubble[1]].add(partial_line)
							else:
								if bubble[0] in seenNode:
									hasDup = True
									#break
								reads_dict[rn_qp].append(bubble[0])
								seenNode.add(bubble[0])
								bubbleCoverByRead[bubble[0]].add(partial_line)
				else:
					if node in seenNode:
						hasDup = True
						break
					reads_dict[rn_qp].append(node)
					seenNode.add(node)
					nodeCoveredByRead[node] += 1
		if hasDup:
			continue	
		for k in range(0,len(reads_dict[rn_qp])-1):
			pair1= str(reads_dict[rn_qp][k])+"_"+str(reads_dict[rn_qp][k+1]) # not taking care of reverse direction now
			pair2= str(reads_dict[rn_qp][k+1])+"_"+str(reads_dict[rn_qp][k])
			# should take of direction, not adding pairs reverse of each other
			if pair1 not in consec_pairs and pair2 not in consec_pairs:
				consec_pairs.add(pair1)
		#if name == 'm150105_192231_42177R_c100761782550000001823161607221526_s1_p0/64358/951_15471':
		#	print('read of interests', reads_dict[rn_qp])

	return reads_dict, consec_pairs, bubbleCoverByRead, nodeCoveredByRead, readnames, chunk_id

def getChunk(gamJson, threadN):
	wc = subprocess.Popen('wc -l %s'%gamJson, shell = True, stdout = subprocess.PIPE)
	stdout, stderr = wc.communicate()
	NR = int(stdout.decode().split()[0])
	if NR % threadN == 0:
		chunkSize = int(NR / threadN)
	else:
		# Why we have a math.py in this folder??? I cannot use math.ceil() now!
		if NR % threadN != 0:
			chunkSize = NR // threadN + 1
		else:
			chunkSize = NR / threadN
		#chunkSize = math.ceil(NR / threadN)
	chunks = []
	for i in range(threadN):
		if (i + 1) * chunkSize <= NR: 
			chunks.append((i*chunkSize, (i+1)*chunkSize))
		else:
			chunks.append((i*chunkSize, NR))
	return chunks, chunkSize

def vg_reader_multithreading(locus_file, gam_file, t):
	reverse_mapping, locus_branch_mapping, allele_reverse_mapping = reverse_map(locus_file)
	for n, b in locus_branch_mapping.items():
		print(n, b)
	bubble_of_interest = -1 #list(reverse_mapping[13153682])[0]
	reads_dict = defaultdict()
	consec_pairs = set()
	bubbleCoverByRead = list()
	nodeCoveredByRead = defaultdict(int)
	readnames = []
	readdict_allele = list()
	for gam in gam_file:
		readnames_i = []
		bubbleCoverByRead_i = defaultdict(set)
		readdict_allele_i = defaultdict(list)
		chunks, chunkSize = getChunk(gam, t)
		print('Reading', gam)
		filecache = open(gam).readlines()
		p = Pool(t)
		results = []
		results_allele = []
		print('Processing', gam)
		for c in range(len(chunks)):
			results.append(p.apply_async(AlignmentParse, args=(filecache[chunks[c][0]:chunks[c][1]], reverse_mapping, c, chunkSize, )))
			results_allele.append(p.apply_async(alleleDetect, args=(filecache[chunks[c][0]:chunks[c][1]], allele_reverse_mapping, c, chunkSize, )))
		p.close()
		p.join()
		for i in results:
			reads_dict_c, consec_pairs_c, bubbleCoverByRead_c_i, nodeCoveredByRead_c, readnames_c_i, c = i.get()
			#print(reads_dict_c)
			readnames_i.append((c, readnames_c_i))
			reads_dict.update(reads_dict_c)
			for pair in consec_pairs_c:
				reverse_pair = pair.split('_')[1]+'_'+pair.split('_')[0]
				if reverse_pair not in consec_pairs:
					consec_pairs.add(pair)
			for b, r in bubbleCoverByRead_c_i.items():
				try:
					bubbleCoverByRead_i[b] = bubbleCoverByRead_i[b].union(r)
				except KeyError:
					bubbleCoverByRead_i[b] = r
			for N, n in nodeCoveredByRead_c.items():
				nodeCoveredByRead[N] += n
		for i in results_allele:
			readdict_allele_c_i = i.get()
			readdict_allele_i.update(readdict_allele_c_i)

		bubbleCoverByRead.append(bubbleCoverByRead_i)
		readnames_i = sorted(readnames_i)
		readnames_i2 = []

		readdict_allele.append(readdict_allele_i)
		for chunk_readnames in readnames_i:
			readnames_i2.extend(chunk_readnames[1])
		readnames.append(readnames_i2)

	return reads_dict, consec_pairs, locus_branch_mapping, bubbleCoverByRead, nodeCoveredByRead, readnames, readdict_allele, bubble_of_interest

#t = 3
#reads_dict, consec_pairs, locus_branch_mapping_raw, bubbleCoverByRead, nodeCoveredByRead = vg_reader_multithreading(sys.argv[1], sys.argv[2:], t)
#print(reads_dict)
