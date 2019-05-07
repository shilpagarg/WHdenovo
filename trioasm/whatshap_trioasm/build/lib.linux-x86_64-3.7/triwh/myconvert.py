#!/usr/bin/python

overlap = 50

def reversecomplement(sequence):
	reverse = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	return ''.join([reverse[nucleotide] for nucleotide in sequence[::-1]])

from collections import defaultdict
import fileinput
import re
import sys
pattern = re.compile(r"^L:([-+]):(\d+):([-+])$")

nodes = []
edges = []

lines = fileinput.input()

number = 0
while True:
	nameline = lines.readline().strip()
	sequence = lines.readline().strip()
	if nameline == "" or sequence == "": break
	nodes.append(sequence)
	parts = nameline.split(' ')
	for part in parts:
		match = pattern.match(part)
		if not match: continue
		numberdirection = match.group(1)
		target = int(match.group(2))
		targetdirection = match.group(3)
		edges.append((number, numberdirection, target, targetdirection))
	number += 1

number = 0
for node in nodes:
	print("S\t" + str(number) + "\t" + node)
	number += 1

for edge in edges:
	leftpart = nodes[edge[0]]
	rightpart = nodes[edge[2]]
	if edge[1] == '-':
		leftpart = reversecomplement(leftpart)
	if edge[3] == '-':
		rightpart = reversecomplement(rightpart)
	leftpart = leftpart[-overlap:]
	rightpart = rightpart[:overlap]
	if leftpart != rightpart:
		sys.stderr.write("incorrect overlap: " + str(edge) + "(" + leftpart + "," + rightpart + ")" + '\n')
	else:
		print('L\t' + str(edge[0]) + '\t' + edge[1] + '\t' + str(edge[2]) + '\t' + edge[3] + '\t' + str(overlap) + 'M')
