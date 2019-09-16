#!/usr/bin/python

from collections import defaultdict
import fileinput
import re

degrees = defaultdict(int)
count = 0
count2 = 0
nodes = set()
edge_nodes = set()
for line in fileinput.input():
	if line[0] == 'S':
		parts2 = line.split('\t')
		nodes.add(parts2[1])
		count += 1
	if line[0] != 'L': 
		continue
	parts = line.split('\t')
	degrees[parts[1]] += 1
	degrees[parts[3]] += 1
	edge_nodes.add(parts[1])
	edge_nodes.add(parts[3])
	count2 += 1
#print(count, count2)

#print(str(sorted([(d, degrees[d]) for d in degrees], key = lambda x: -x[1])).replace('), ', ')\n').translate("'[](),\''", ""))
print(str(sorted([(d, degrees[d]) for d in degrees], key = lambda x: -x[1])).replace('), ', ')\n').translate(None, '[](),\''))
#print(str(sorted([(d, degrees[d]) for d in degrees], key = lambda x: -x[1])).replace('), ', ')\n'))
#nodes.append('10000000')
#print(nodes)
for element in nodes:
	if element not in edge_nodes:
		print((element + " 0").replace('), ', ')').translate(None, '[](),\''))