#!/usr/bin/python

import fileinput

currentread = ""
nodeid = 1
node_size = 500

def output_read(seq):
	seqs = seq.split("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
	for s in seqs:
		output_one_read(s)

def read_to_nodegroups(seq):
	from cStringIO import StringIO
	current_seq = StringIO()
	current_len = 0
	for i in range(0, len(seq)):
		if seq[i] in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']:
			current_seq.write(seq[i])
			current_len += 1
			if current_len > node_size:
				yield [current_seq.getvalue()]
				current_seq = StringIO()
				current_len = 0
		else:
			yield [current_seq.getvalue()]
			current_seq = StringIO()
			if seq[i] == 'N':
				yield ['A', 'C', 'G', 'T']
			elif seq[i] == 'R':
				yield ['A', 'G']
			elif seq[i] == 'Y':
				yield ['A', 'C', 'T']
			elif seq[i] == 'K':
				yield ['A', 'G', 'T']
			elif seq[i] == 'M':
				yield ['A', 'C']
			elif seq[i] == 'S':
				yield ['C', 'G']
			elif seq[i] == 'W':
				yield ['A', 'T']
			elif seq[i] == 'B':
				yield ['C', 'G', 'T']
			elif seq[i] == 'D':
				yield ['A', 'G', 'T']
			elif seq[i] == 'H':
				yield ['A', 'C', 'T']
			elif seq[i] == 'V':
				yield ['A', 'C', 'G']
	last = current_seq.getvalue()
	if len(last) > 0: yield [last]

def filter_start_N(seq):
	for i in range(0, len(seq)):
		if seq[i] != 'N':
			seq = seq[i:]
			break
	for i in range(len(seq)-1, 0, -1):
		if seq[i] != 'N':
			seq = seq[:i]
			break
	return seq

def output_groups(seq):
	global nodeid
	if len(seq) == 0: return
	last_nodes = []
	for group in read_to_nodegroups(seq):
		new_nodes = []
		for seq in group:
			print("S\t" + str(nodeid) + "\t" + seq)
			new_nodes.append(nodeid)
			for oldnode in last_nodes:
				print("L\t" + str(oldnode) + "\t+\t" + str(nodeid) + "\t+\t0M")
			nodeid += 1
		last_nodes = new_nodes

def output_one_read(seq):
	seq = filter_start_N(seq)
	output_groups(seq)


for line in fileinput.input():
	l = line.strip()
	if l[0] == '>':
		output_read(currentread)
		currentread = ""
	else:
		currentread += l

