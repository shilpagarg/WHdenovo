import sys
#from Bio import SeqIO

#all_partitioned_readnames all_reads.fq
partitioned_reads = set()
total_reads = set()
with open(sys.argv[1], 'r') as f:
	for line in f:
		var = line.rstrip()
		partitioned_reads.add(str(var))


with open(sys.argv[2], 'r') as f:
	for line in f:
		if line[0] == ">":
			total_reads.add(line[1:].rstrip())
#print(total_reads)
unpartitioned_reads = total_reads - partitioned_reads
print(len(unpartitioned_reads))

with open(sys.argv[1]+'.unpartitioned', 'w') as f:
	for element in unpartitioned_reads:
		f.write(element + "\n")
