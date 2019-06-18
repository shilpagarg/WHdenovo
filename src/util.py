from collections import defaultdict

class aux_unitigs():
	# a module for filtering the unitig connections
	def __init__(self, unitigs):
		self.unitigs = unitigs
		self.endMap = dict()
		self.unitigStart = dict()
		self.unitigEnd = dict()
		for u in range(len(unitigs)):
			self.endMap[unitigs[u][-1]] = u
			self.endMap[unitigs[u][0]] = u
			self.unitigStart[unitigs[u][-1]] = unitigs[u][0]
			self.unitigEnd[unitigs[u][0]] = unitigs[u][-1]
		self.connections = defaultdict(int)
		self.pair = defaultdict(set)
	def addConnection(self, new_g):
		# add connections detected from reads. Merge equivalent 
		# connections together
		connections = []
		for i in range(0, len(new_g), 2):
			node = new_g[i]
			if node in self.unitigStart:
				connections.append((self.endMap[node], -1))
			if node in self.unitigEnd:
				connections.append((self.endMap[node], 1))
		for i in range(len(connections) - 1):
			u1 = connections[i]
			u2 = connections[i + 1]
			if (u1[0], u2[0]) in self.pair:
				if (u1, u2) in self.connections:
					self.connections[(u1, u2)] += 1
				elif equivConn((u1, u2)) in self.connections:
					self.connections[equivConn((u1, u2))] += 1
				else:
					self.connections[(u1, u2)] = 1
					self.pair[u1[0], u2[0]].add((u1, u2))
			elif (u2[0], u1[0]) in self.pair:
				if (u1, u2) in self.connections:
					self.connections[(u1, u2)] += 1
				elif equivConn((u1, u2)) in self.connections:
					self.connections[equivConn((u1, u2))] += 1
				else:
					self.connections[(u1, u2)] = 1
					self.pair[(u2[0], u1[0])].add((u1, u2))
			else:
				self.pair[(u1[0], u2[0])].add((u1, u2))
				self.connections[(u1, u2)] += 1

	def returnConfident(self):
		# First pick the most prevalent type of connection between 
		# two specific unitigs. Then sort all connections by number 
		# of supporting reads and pick from top, allowing only one 
		# connection for one unitig end.
		allout = []
		for pair, connections in self.pair.items():
			out = (0, ())
			for connection in connections:
				if self.connections[connection] > out[0]:
					out = (self.connections[connection], connection)
			allout.append(out)
		allout = sorted(allout)[::-1]
		AUX = []
		connected = set()
		for connection in allout:
			# connection = (n, (uid1, dir1), (uid2, dir2))
			# n is the frequency of this type of connection. 
			new_g = []
			if connection[0] < 5:
				continue
			SE1 = (connection[1][0][0], connection[1][0][1], 0)
			eSE1 = equivStickyEnd(SE1)
			SE2 = (connection[1][1][0], connection[1][1][1], 1)
			eSE2 = equivStickyEnd(SE2)
			if SE1 in connected or eSE1 in connected or SE2 in connected or eSE2 in connected:
				#if connected[connection[1][0][0]] == 2 or connected[connection[1][1][0]] == 2:
				continue
			print(connection)
			connected.add(SE1)
			connected.add(SE2)
			#connected[connection[1][0][0]] += 1
			#connected[connection[1][1][0]] += 1
			if connection[1][0][1] == 1:
				new_g.append(self.unitigs[connection[1][0][0]][0])
				new_g.append(self.unitigs[connection[1][0][0]][-1])
			else:
				new_g.append(self.unitigs[connection[1][0][0]][-1])
				new_g.append(self.unitigs[connection[1][0][0]][0])
			if connection[1][1][1] == 1:
				new_g.append(self.unitigs[connection[1][1][0]][0])
				new_g.append(self.unitigs[connection[1][1][0]][-1])
			else:
				new_g.append(self.unitigs[connection[1][1][0]][-1])
				new_g.append(self.unitigs[connection[1][1][0]][0])
			AUX.append(new_g)
		return AUX

def equivConn(connection):
	u1 = connection[0]
	u2 = connection[1]
	return ((u2[0], -1 * u2[1]), (u1[0], -1 * u1[1]))

def equivStickyEnd(end):
	# if in a connection we see ((0, 1), (1, 1))
	# Then I define "stickyEnd" for this pair as (0, 1, 0) and (1, 1, 1), 
	# and their equivalents are (0, -1, 1) and (1, -1, 0), 
	# where the 1st of the triple is unitig ID, 2nd is the direction against 
	# the unitig's original direction (1 for same -1 for reverse), and in the
	# 3rd position 0 for on the left 1 for on the right.
	return (end[0], -1 * end[1], abs(end[2] - 1))
