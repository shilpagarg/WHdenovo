import sys
import stream
import logging
import vg_pb2 
from collections import Counter
from collections import defaultdict
import collections
from collections import OrderedDict
from itertools import groupby

locus_file = sys.argv[1]
contigs_input = sys.argv[2]
#trans = sys.argv[3]


count=0

locus_count=0
prev_startsnarl = 0
prev_endsnarl = 0
locus_branch_mapping=OrderedDict()
locus_count=0
prev_startsnarl = 0
prev_startsnarl_orientation = -1
prev_endsnarl = 0
prev_endsnarl_orientation = -1
reads_dict = defaultdict(list)
path_in_bubble =[]
with stream.open(str(locus_file), "rb") as istream:
	for data in istream:
		l = vg_pb2.SnarlTraversal()
		l.ParseFromString(data)
		#TODO: make ordered doctionary locus_branch_mapping
		# handle forward and backward case of nodes
		current_startsnarl = l.snarl.start.node_id
		current_startsnarl_orientation = l.snarl.start.backward
		current_endsnarl = l.snarl.end.node_id
		current_endsnarl_orientation = l.snarl.end.backward

		if current_startsnarl == prev_startsnarl and current_endsnarl == prev_endsnarl:
			path_in_bubble.append(l)
		else:
			
			locus_branch_mapping[locus_count]=path_in_bubble
			locus_count=locus_count-1
			#per_locus = []
			#per_locus.append(path_in_bubble)
			path_in_bubble = []
			path_in_bubble.append(l)
			
		prev_startsnarl = current_startsnarl
		prev_startsnarl_orientation = current_startsnarl_orientation
		prev_endsnarl = current_endsnarl
		prev_endsnarl_orientation = current_endsnarl_orientation

		#locus_branch_mapping[locus_count]=per_locus[0]

#print(locus_branch_mapping)
#ostream = stream.open(str(trans), 'wb')
print(locus_branch_mapping)
count=0
name = str(sys.argv[2]).replace('.final_ctgs', '')
with open(contigs_input) as fp:
	for line in fp:
		count+=1
		tokens = line.rstrip().split(",")
		if len(tokens) > 2:
			ostream = stream.open(str(name)+'_contigs_' +str(count)+'.trans', 'wb')
		#for i in range(0, len(tokens)-1):
		#	if int(tokens[i]) in locus_branch_mapping:
		#		ostream = stream.open('contigs_'+str(count)+'.trans', 'wb')
			for i in range(0, len(tokens)-1):
				if int(tokens[i]) in locus_branch_mapping:
					for k in locus_branch_mapping[int(tokens[i])]:
					#	print("HELLO", int(tokens[i]))
						ostream.write(k)

