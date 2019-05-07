import networkx
import sys   
from collections import defaultdict
import stream

inp = sys.argv[1]
out = sys.argv[2] 
outf = open(out, 'w')
   
#This class represents a directed graph  
# using adjacency list representation 
class Graph: 
   
	def __init__(self,vertices): 
		#No. of vertices 
		self.V= vertices  
		  
		# default dictionary to store graph 
		self.graph = defaultdict(list)  
   
	# function to add an edge to graph 
	def addEdge(self,u,v): 
		self.graph[u].append(v) 
   
	'''A recursive function to print all paths from 'u' to 'd'. 
	visited[] keeps track of vertices in current path. 
	path[] stores actual vertices and path_index is current 
	index in path[]'''
	def printAllPathsUtil(self, u, d, visited, path): 
  
		# Mark the current node as visited and store in path 
		visited[u]= True
		path.append(u) 
  
		# If current vertex is same as destination, then print 
		# current path[] 
		if u ==d: 
			print(path) 
		else: 
			# If current vertex is not destination 
			#Recur for all the vertices adjacent to this vertex 
			for i in self.graph[u]: 
				if visited[i]==False: 
					self.printAllPathsUtil(i, d, visited, path) 
					  
		# Remove current vertex from path[] and mark it as unvisited 
		path.pop() 
		visited[u]= False
   
   
	# Prints all paths from 's' to 'd' 
	def printAllPaths(self,s, d): 
  
		# Mark all the vertices as not visited 
		visited =[False]*(self.V) 
  
		# Create an array to store paths 
		path = [] 
  
		# Call the recursive helper function to print all paths 
		path = self.printAllPathsUtil(s, d,visited, path) 
		return path


def vg_graph_reader(inp):
	with stream.open(str(inp), "rb") as istream:
		for data in istream:
			l = vg_pb2.Graph()
			l.ParseFromString(data)
			g = Graph(len(l.node))
			for j in range(len(l.edge)):
				from_edge = getattr(l.edge[j], "from")
				g.addEdge(from_edge, l.edge[j].to)

	return g

g = vg_graph_reader(inp)
num_of_vertices = len(g.V)


for i in range(0,20):
	ran1 = random.randint(0, num_of_vertices)
	ran2 = random.randint(0, num_of_vertices)
	paths = g.printAllPaths(g.V[ran1], g.V[ran2])
	for k in paths:
		outf.write(k)

