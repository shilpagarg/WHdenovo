import sys
import stream
import vg_pb2
from collections import defaultdict
import time
file_input = sys.argv[1]
#start = time.clock()

def reverse_complement(seq):
        seq_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'T', 'c': 'G', 'g': 'C', 't': 'A'}
        return "".join([seq_dict[base] for base in reversed(seq)])

#nodes_list = set()

orderalignment = defaultdict(list)
#orderalignment = defaultdict(lambda: [-1]*100, orderalignment)
with stream.open(str(file_input), "rb") as istream:
        for data in istream:
                g = vg_pb2.Alignment()
                g.ParseFromString(data)
                canu_name = g.name
                canu_chunk_num = int(g.query_position)
                orderalignment[canu_name].append([canu_chunk_num, [g]])
#for i in orderalignment["S1_43_mom1"]:
#        print(i[0], i[1][0].name)
#exit()
#new_orderalignment = defaultdict(list)
for k,v in orderalignment.items():
        if len(v) > 1:
                #print(k, orderalignment[k])
                v.sort()
                orderalignment[k] = v

                #print(k, orderalignment[k])
        #new_orderalignment[k] = [x for x in v if x != -1]

ostream = stream.open(sys.argv[2], 'wb')

for k,v in orderalignment.items():
        new_g = vg_pb2.Alignment()
        new_g.name = v[0][1][0].name
        for i in range(len(v)):
                #new_g.name = v[i][1][0].name
                #tmp = v[i].name.split("_")[0]
                for j in range(len(v[i][1][0].path.mapping)):
                        new_g.path.mapping.extend([v[i][1][0].path.mapping[j]])
#ostream1 = stream.open(file_input + '_combined.gam', 'wb')
#ostream1.write(new_g)
#ostream1.close()
        ostream.write(new_g)


#end = time.clock()
#print('TIME:', end - start)
