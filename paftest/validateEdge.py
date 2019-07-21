import sys

"""
Usage: python validateEdge.py GFA truth1 truth2 > wrongEdges.gfa
"""

gfa = open(sys.argv[1], 'r').read().split('\n')[:-1]
p1 = set(open(sys.argv[2], 'r').read().split('\n')[:-1])
p2 = set(open(sys.argv[3], 'r').read().split('\n')[:-1])


total_edge = 0
correct = 0
wrong = 0
unknown = 0
for line in gfa:
	if line[0] != 'L':
		continue
	tokens = line.split('\t')
	if tokens[-1] == 'WRONG' or tokens[-1] == 'DASH' or tokens[-1] == 'ONETWO':
		continue
	total_edge += 1
	r1 = tokens[1]
	r2 = tokens[3]
	if r1 in p1 and r2 in p2:
		if tokens[-1] == 'RED':
			correct += 1
		else:
			wrong += 1
			print(line)
	elif r1 in p2 and r2 in p1:
		if tokens[-1] == 'RED':
			correct += 1
		else:
			wrong += 1
			print(line)
	elif r1 in p1 and r2 in p1:
		if tokens[-1] == 'GREEN':
			correct += 1
		else:
			wrong += 1
			print(line)
	elif r1 in p2 and r2 in p2:
		if tokens[-1] == 'GREEN':
			correct += 1
		else:
			wrong += 1
			print(line)
	else:
		unknown += 1

sys.stderr.write('Correct edge:' + str(correct)+ '\n')
sys.stderr.write('Wrong edge:' + str(wrong) + '\n')
sys.stderr.write('Unknown edge:' + str(unknown) + '\n')
sys.stderr.write('Total edge:' + str(total_edge)+ '\n')
sys.stderr.write('correct/known:'+str(correct/(total_edge-unknown))+ '\n')
