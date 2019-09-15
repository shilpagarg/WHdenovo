from collections import defaultdict, OrderedDict
from multiprocessing import Pool, Process, Manager
from operator import attrgetter
import subprocess
import sys
import time
import math
import json
import stream
import vg_pb2

def ceil(n, pos = True):
	# whatshap has its own "math" package in the src directory.
	# I cannot use the python original math library.
	# pos is for when you want the result to be positive (at least 1)
	# ceil(0) -> 1, ceil(1.1) > 2, ceil(-1) -> 1, ceil(-1.1, pos = False) -> -1
	if n % 1 == 0:
		if n > 0:
			return int(n)
		else:
			if not pos:
				return int(n)
			else:
				return 1
	else:
		if n // 1 + 1 < 1:
			if not pos:
				return int(n // 1 + 1)
			else:
				return 1
		else:
			return int(n // 1 + 1)

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
	prev_startsnarl_orientation = -1
	prev_endsnarl_orientation = -1
	insidebubble = 0
	with stream.open(str(locus_file), "rb") as istream:
		for data in istream:
			l = vg_pb2.SnarlTraversal()
			l.ParseFromString(data)
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
					pass
				else:
					try:
						locus_branch_mapping[locus_count] = per_locus
					except NameError:
						pass
					locus_count -= 1
					per_locus = []
			else:
				if current_startsnarl == prev_startsnarl and current_endsnarl == prev_endsnarl and current_endsnarl_orientation == prev_endsnarl_orientation and prev_startsnarl_orientation == current_startsnarl_orientation:
					if insidebubble == 2:
						path_in_bubble = mergePath(tempPath, path_in_bubble, insideBack, pathBack, local_path_back)
						per_locus.append(path_in_bubble)
						insidebubble = 0
						insideBack = False
						pathBack = False
					else:
						per_locus.append(path_in_bubble)
				else:
					if insidebubble == 1:
						insidebubble = 2
						path_in_bubble = mergePath(tempPath, path_in_bubble, insideBack, pathBack, local_path_back)
						per_locus.append(path_in_bubble)
					else:
						try:
							locus_branch_mapping[locus_count] = per_locus
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
	het_count= 0
	alleles_per_pos = dict()
	for k,bubble in locus_branch_mapping.items():
		alleles_per_pos[k] = len(bubble)
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
	return reverse_mapping, allele_reverse_mapping, alleles_per_pos, locus_branch_mapping

class alignmentSets(object):
	# An object for storing partial alignments from separate individuals.
	# It does the post-processing on them, including global alignment, repeat
	# detection and removing, edge building
	def __init__(self, mode = 'trio'):
		if mode == 'trio':
			self.partialList = [[], [], []]
			self.fullReadList = [[], [], []]
			self.bubbleReadMap = [defaultdict(set), defaultdict(set), defaultdict(set)]
		elif mode == 'individual':
			self.partialList = [[]]
			self.fullReadList = [[]]
			self.bubbleReadMap = [defaultdict(set)]

	def addPartial(self, partial, sample = 0):
		self.partialList[sample].append(partial)

	def mergeChunk(self, alignmentSetsObj):
		for i in range(len(self.partialList)):
			self.partialList[i].extend(alignmentSetsObj.partialList[i])

	def postProcessing(self, lowc, highc):
		for sample in range(len(self.partialList)):
			self.globalAlign(sample)
		self.getBubbleReadMap()
		self.trimRepeat(lowc, highc)
		self.getBubbleReadMap()

	def globalAlign(self, sample):
		# Gather all the partials of one read name together in the order of 
		# query_position.
		self.partialList[sample] = sorted(self.partialList[sample], key=attrgetter('name', 'query_position'))
		lastname = ''
		for partial in self.partialList[sample]:
			if partial.name != lastname:
				if lastname != '':
					self.fullReadList[sample].append(read)
				read = fullRead()
				read.addPartial(partial)
			else:
				
				read.addPartial(partial)
			lastname = partial.name
		self.fullReadList[sample].append(read)
		print('Sample', sample, 'total read number', len(self.fullReadList[sample]))
					
	def depth(self, node):
		cov = 0
		for sample in range(len(self.bubbleReadMap)):
			cov += len(self.bubbleReadMap[sample][node])
		return cov

	def removeNodeFromRead(self, node, read):
		for partial in read.partials:
			pL = len(partial)
			for i in range(pL):
				if partial[i] == node:
					partial[i] = None
			for i in range(pL):
				try:
					partial.remove(None)
				except:
					break
		if node < 0:
			aL = len(read.alleles)
			for i in range(aL):
				if read.alleles[i][0] == node:
					read.alleles[i] = None
			for i in range(aL):
				try:
					read.alleles.remove(None)
				except:
					break

	def trimRepeat(self, lowc, highc):
		print('Start detecting repeated nodes/bubbles')
		repeatCollection = defaultdict(int)
		repeatToRead = defaultdict(list)
		for sample in range(len(self.fullReadList)):
			for read in self.fullReadList[sample]:
				seenNode = defaultdict(int)
				for partial in read.partials:
					for node in partial:
						seenNode[node] += 1
				for n, c in seenNode.items():
					if c >= 2:
						repeatCollection[n] += 1
						repeatToRead[n].append((read, sample))
		print('Repeat information collected')
		print('Totally %d repetitive nodes/bubbles found'%len(list(repeatCollection.keys())))
		count = 0
		child_count = 0
		threshold = 0.3
		print('lowc:', lowc, 'highc:', highc)
		for n, c in repeatCollection.items():
			cov = self.depth(n)
			if c / cov < threshold and cov <= highc and cov >= lowc:
				for (read, sample) in repeatToRead[n]:
					try:
						self.fullReadList[sample].remove(read)
						print('Removing read', read)
						count += 1
						if sample == 2:
							child_count += 1
					except ValueError:
						pass
					
		print('Totally %d reads removed' % count)
		print('Totally %d child reads removed' % child_count)
		del repeatToRead
		for sample in range(len(self.fullReadList)):
			for read in self.fullReadList[sample]:
				for partial in read.partials:
					pL = len(partial)
					for i in range(pL):
						cov = self.depth(partial[i])
						if repeatCollection[partial[i]] / cov >= threshold or cov <= lowc or cov >= highc:
							partial[i] = None
					for i in range(pL):
						try:
							partial.remove(None)
						except:
							break
				aL = len(read.alleles)
				for var in range(aL):
					cov = self.depth(read.alleles[var][0])
					if repeatCollection[read.alleles[var][0]] / cov >= threshold or cov <= lowc or cov >= highc:
						read.alleles[var] = None
				for var in range(aL):
					try:
						read.alleles.remove(None)
					except:
						break
	def getBubbleReadMap(self):
		if len(self.bubbleReadMap) == 3:
			self.bubbleReadMap = [defaultdict(set), defaultdict(set), defaultdict(set)]
		if len(self.bubbleReadMap) == 1:
			self.bubbleReadMap = [defaultdict(set)]
		for sample in range(len(self.fullReadList)):
			bubblereadlist = set()
			for read in self.fullReadList[sample]:
				for partial in read.partials:
					for node in partial:
						self.bubbleReadMap[sample][node].add(read)
						if node < 0:
							bubblereadlist.add(read)
			print('Sample', sample, 'has', len(bubblereadlist), 'reads associated with bubble')

	def getEdges(self):
		edges = set()
		for sample in range(len(self.fullReadList)):
			for read in self.fullReadList[sample]:
				for partial in read.partials:
					if len(partial) > 1:
						for i in range(len(partial)-1):
							edge1 = (partial[i], partial[i + 1])
							edge2 = (partial[i + 1], partial[i])
							if edge1 not in edges:
								edges.add(edge2)
		return edges

def coverage(bRM, node):
	cov = 0
	for sample in range(len(bRM)):
		cov += bRM[sample][node]
	return cov

class fullRead(object):
	# Object for storing read information after global alignment
	def __init__(self):
		self.name = None
		self.partials = []
		self.alleles = []
		self.query_pos = []
	def addPartial(self, partial):
		if partial.query_position in self.query_pos:
			pass
		else:
			if self.name != None:
				assert partial.name == self.name
			self.name = partial.name
			self.partials.append(partial.path)
			self.alleles.extend(partial.alleles)		
	def __str__(self):
		return self.name
	def __repr__(self):
		return self.name
class partial(object):
	# Simply parse a line of json content, and store the information.
	def __init__(self, jsonLine, reverse_mapping, allele_reverse_mapping):
		self.name = None
		self.score = 0
		self.query_position = 0
		self.rawMapping = []
		self.path = []
		self.alleles = []
		self.parseJson(jsonLine)
		self.getPath(reverse_mapping)
		self.getAllele(allele_reverse_mapping)
		del self.rawMapping
	def __str__(self):
		return self.name
	def __repr__(self):
		return self.name
	def parseJson(self, jsonLine):
		g = json.loads(jsonLine)
		self.name = g['name']
		try:
			self.query_position = int(g['query_position'])
		except KeyError:
			self.query_position = 0
		try:
			self.score = int(g['score'])
		except KeyError:
			self.score = 0
		self.rawMapping = g['path']['mapping']

	def getPath(self, reverse_mapping):
		# Get a node based path.
		if len(self.rawMapping) == 1:
			node = int(self.rawMapping[0]['position']['node_id'])
			if node in reverse_mapping:
				bubble = list(reverse_mapping[node])
				self.path = bubble
			else:
				self.path.append(node)
		else:
			for i in range(len(self.rawMapping)):
				node = int(self.rawMapping[i]['position']['node_id'])
				if node in reverse_mapping:
					bubble = reverse_mapping[node]
					if len(bubble) == 1:
						try:
							if self.path[-1] == list(bubble)[0]:
								continue
						except IndexError:
							pass
						self.path.append(list(bubble)[0])
					else:
						bubble = list(bubble)
						if len(self.path) == 0:
							nextnode = int(self.rawMapping[i+1]['position']['node_id'])
							self.path.append(bubble[0])
							self.path.append(bubble[1])
							if list(reverse_mapping[nextnode])[0] == bubble[0]:
								self.path.reverse()
						else:
							if self.path[-1] == bubble[0]:
								self.path.append(bubble[1])
							else:
								self.path.append(bubble[0])
				else:
					self.path.append(node)

	def getAllele(self, allele_reverse_mapping):
		# Get an allele based variant sequence
		prev_tmp = []
		prev_locus = -1
		added = False
		for i in range(len(self.rawMapping) - 1):
			edge1 = (int(self.rawMapping[i]['position']['node_id']), int(self.rawMapping[i+1]['position']['node_id']))
			edge2 = (int(self.rawMapping[i+1]['position']['node_id']), int(self.rawMapping[i]['position']['node_id']))
			if edge1 in allele_reverse_mapping or edge2 in allele_reverse_mapping:
				if edge1 in allele_reverse_mapping:   
					#qualities = [10]* reverse_mapping[edge1][0][2]
					qualitie = 1 
					node_inf = [tuple(i[0:3]) for i in allele_reverse_mapping[edge1]]
				else:
					# qualities = [10]* reverse_mapping[edge2][0][2]
					qualities = 1 
					node_inf = [tuple(i[0:3]) for i in allele_reverse_mapping[edge2]]
				tmp = node_inf.copy()
				if prev_locus != tmp[0][0]:
					added = False
					prev_tmp = tmp.copy()
					prev_locus = tmp[0][0]
				if added:
					continue
				interset_tmp = list(set(tmp).intersection(set(prev_tmp)))
				if len(interset_tmp) == 1:
					qualities = 1 
					self.alleles.append((interset_tmp[0][0], interset_tmp[0][1]))
					added = True

def runChunk(jsonLines, sample, reverse_mapping, allele_reverse_mapping, mode, results):
	
	alnSet = alignmentSets(mode = mode)
	for line in jsonLines:
		alnSet.addPartial(partial(line, reverse_mapping, allele_reverse_mapping), sample)
	results.append(alnSet)

def getChunk(gamFilecache, threadN):
	NR = len(gamFilecache)
	if NR % threadN == 0:
		chunkSize = 10000
	else:
		if NR % threadN != 0:
			chunkSize = NR // threadN + 1
		else:
			chunkSize = NR / threadN
	chunkSize = 10000
	chunks = []
	for i in range(threadN):
		if (i + 1) * chunkSize <= NR: 
			chunks.append((i*chunkSize, (i+1)*chunkSize))
		else:
			chunks.append((i*chunkSize, NR))
	return chunks

def vg_read(locus_file, gam_file, t, lowc, highc):
	reverse_mapping, allele_reverse_mapping, alleles_per_pos, locus_branch_mapping = reverse_map(locus_file)
	if len(gam_file) == 3:
		mode = 'trio'
	elif len(gam_file) == 1:
		mode = 'individual'
	totalAlnSet = alignmentSets(mode = mode)
	for i in range(len(gam_file)):
		print('Reading', gam_file[i])
		filecache = open(gam_file[i]).readlines()
		chunkidx = getChunk(filecache, t)
		m = Manager()
		results = m.list()
		p = Pool(t)
		processes = []
		for c in range(len(chunkidx)):
			processes.append(p.apply_async(runChunk, (filecache[chunkidx[c][0]:chunkidx[c][1]], i, reverse_mapping, allele_reverse_mapping, mode, results)))
			#p = Process(target = runChunk, args = (filecache[chunkidx[c][0]:chunkidx[c][1]], i, reverse_mapping, allele_reverse_mapping, mode, results))
			#p.start()
			#pool.append(p)
		#for p in pool:
		#	p.join()
		p.close()
		p.join()
		for i in results:
			totalAlnSet.mergeChunk(i)
	totalAlnSet.postProcessing(lowc, highc)
	edges = totalAlnSet.getEdges()

	return totalAlnSet, edges, alleles_per_pos, locus_branch_mapping
