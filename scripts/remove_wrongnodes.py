import sys

wrongnodes = set()
with open(sys.argv[1], 'r') as m:
	for line in m:
		var = line.rstrip()
		wrongnodes.add(var)
m.close()

g = open(sys.argv[3], 'w')
print(wrongnodes)
with open(sys.argv[2], 'r') as f:
	for line in f:
		var = line.rstrip().split("\t")
		if var[0] == "S":
			if var[1] in wrongnodes:
				#print('yay')
				continue
		elif var[0] == "L":
			if (var[1] in wrongnodes) or (var[3] in wrongnodes):
				#print('yay2')
				continue
		g.write(line)
