from collections import defaultdict

class aux_unitigs():
	"""docstring for ClassName"""
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
		allout = []
		for pair, connections in self.pair.items():
			out = (0, ())
			for connection in connections:
				if self.connections[connection] > out[0]:
					out = (self.connections[connection], connection)
			#if out[1] >= 100:
			allout.append(out)
		allout = sorted(allout)[::-1]
		AUX = []
		connected = defaultdict(int)
		for connection in allout:
			new_g = []
			if connected[connection[1][0][0]] == 2 or connected[connection[1][1][0]] == 2:
				continue
			connected[connection[1][0][0]] += 1
			connected[connection[1][1][0]] += 1
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
'''
unitigs = [[1,2,3],[4,5,6],[7,8,9]]

aux = aux_unitigs(unitigs)
aux.addConnection([3,1,4,6,9,7]) # 0,-1  1,1  2,-1
aux.addConnection([6,4,1,3])     # 1,-1  0,1
aux.addConnection([1,3, 4, 6])   # 0,1   1,1
aux.addConnection([7, 9, 6,4, 1, 3])   # 0,1   1,1
aux.addConnection([6,4,3,1])     # 1,-1  0,-1
aux.addConnection([1,3,6,4])
aux.addConnection([4,6,3,1])
aux.addConnection([6,4,9,7])
print('==================')
print(aux.connections)
print(aux.pair)
AUX = aux.returnConfident()
print(AUX)
'''
