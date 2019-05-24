import sys

graph = open(sys.argv[1], 'r').read().split('\n')[:-1]

total_edge = 0
error_edge = 0
dash_edge = 0
wrong_edge = 0
node_count = 0
partitioned = set()
wrong = open('wrongColor.txt','w')
for line in graph:
	if line[0] == 'S':
		node_count += 1
		continue
	total_edge += 1
	tokens = line.split('\t')
	read1 = tokens[1]
	read2 = tokens[3]
	color = tokens[6]
	if color == 'GREEN':
		partitioned.add(read1)
		partitioned.add(read2)
	#if read1[-1] == read2[-1] and color == 'RED':
	#	wrong.write(line+'\n')
	#	error_edge += 1
	if read1[-1] != read2[-1] and color == 'GREEN':
		wrong.write(line+'\n')
		error_edge += 1
	if color == 'DASH':
		dash_edge += 1
	if color == 'somethingWrong':
		wrong_edge += 1

print('Total edges: %d'%total_edge)
print('Total nodes: %d'%node_count)
print('Nodes connected by GREEN: %d'%len(partitioned))
print('GREEN Error rate: %f' % (error_edge/total_edge))
print('DASH edge count: %d'%dash_edge)
print('BUG edge: %d'%wrong_edge)
