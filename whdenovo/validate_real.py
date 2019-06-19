from collections import defaultdict
import sys, os, subprocess


def readBothFile(allreads, tagTXT):
	partitioned = open(allreads, 'r').read().split('\n')[:-1]
	readnameSet1 = set()
	for read in partitioned:
		readnameSet1.add(read.split()[0])
	
	tag = open(tagTXT, 'r').read().split('\n')[:-2]
	readnameSet2 = set()
	tagged = []
	for read in tag:
		info = read.split('\t')
		if info[1] != '':
			readnameSet2.add(info[0])
			tagged.append(read)
	
	readIsec = readnameSet1.intersection(readnameSet2)

	return partitioned, tagged, readIsec

def partition_selection(partitioned, readIsec):
	read2block = defaultdict(list)
	for line in partitioned:
		readname = line.split(' ')[0]
		nVar = int(line.split(' ')[3])
		if readname in readIsec:
			read2block[readname].append((nVar, line))
	selected = []
	for k, v in read2block.items():
		v = sorted(v)
		selected.append(v[-1][1])
	return selected

def readPred(partitioned, readIsec):
	pred_hap1 = defaultdict(list)
	pred_hap2 = defaultdict(list)
	considered_total = 0
	# Pick the record with more variants.
	read_block_hp_nvar = defaultdict(list)
	for read in partitioned:
		tokens = read.split()
		if tokens[0] in readIsec:
			read_block_hp_nvar[tokens[0]].append((tokens[1], tokens[2], int(tokens[3])))
	for rn, block_hp_nvar in read_block_hp_nvar.items():
		# Pick the record with more variants.
		block_hp_nvar = sorted(block_hp_nvar, key = lambda t:t[2])[-1] # returning a tuple here
		if block_hp_nvar[1] == '1':
			pred_hap1[block_hp_nvar[0]].append(rn)
		else:
			pred_hap2[block_hp_nvar[0]].append(rn)
		considered_total += 1
	return pred_hap1, pred_hap2, considered_total, read_block_hp_nvar

def readTrue(tagged, readIsec):
	trueHap1 = []#defaultdict(list)
	trueHap2 = []#defaultdict(list)
	
	for read in tagged:
		info = read.split()
		readname = info[0]
		if readname in readIsec:
			trueHap = int(info[1].split(':')[-1])
			blockID = info[2]
			if trueHap == 1:
				trueHap1.append(readname)
				#trueHap1[blockID].append(readname)
			elif trueHap == 2:
				trueHap2.append(readname)
				#trueHap2[blockID].append(readname)

	return trueHap1, trueHap2 # For totally how many reads we have


def compute_accuracy(pred_hap1, pred_hap2, trueHap1, trueHap2, read_block_hp_nvar):
	predBlocks = list(pred_hap1.keys())
	successCount = 0
	correct_reads = []
	wrong_reads = []
	for k1 in predBlocks:
			successCount += max(len(list(set(pred_hap1[k1]).intersection(set(trueHap1)))), len(list(set(pred_hap1[k1]).intersection(set(trueHap2)))))
			successCount += max(len(list(set(pred_hap2[k1]).intersection(set(trueHap1)))), len(list(set(pred_hap2[k1]).intersection(set(trueHap2)))))
			if len(list(set(pred_hap1[k1]).intersection(set(trueHap1)))) > len(list(set(pred_hap1[k1]).intersection(set(trueHap2)))):
				correct_reads.extend(list(set(pred_hap1[k1]).intersection(set(trueHap1))))
				wrong_reads.extend(list(set(pred_hap1[k1]).difference(set(trueHap1))))
			else:
				correct_reads.extend(list(set(pred_hap1[k1]).intersection(set(trueHap2))))
				wrong_reads.extend(list(set(pred_hap1[k1]).difference(set(trueHap2))))
			if len(list(set(pred_hap2[k1]).intersection(set(trueHap1)))) > len(list(set(pred_hap2[k1]).intersection(set(trueHap2)))):
				correct_reads.extend(list(set(pred_hap2[k1]).intersection(set(trueHap1))))
				wrong_reads.extend(list(set(pred_hap2[k1]).difference(set(trueHap1))))
			else:
				correct_reads.extend(list(set(pred_hap2[k1]).intersection(set(trueHap2))))
				wrong_reads.extend(list(set(pred_hap2[k1]).difference(set(trueHap2))))

	partitioned_accuracy = successCount / considered_total
	success_rate = successCount / len(tagged)
	print('partitioned:', considered_total)
	print('tagged:', len(tagged))
	print('intersection:', len(readIsec))
	#print('partitioned/tagged:', considered_total/len(tagged))
	print('successCount/consideredPartitioned:', partitioned_accuracy)
	oc = open('correct.reads', 'w')
	for i in correct_reads:
		bhn = read_block_hp_nvar[i]
		bhn = sorted(bhn, key = lambda t:t[2])[-1]
		hp = bhn[1]
		oc.write(i + '\t' + hp + '\n')
	oc.close()
	ow = open('wrong.reads', 'w')
	for i in wrong_reads:
		bhn = read_block_hp_nvar[i]
		bhn = sorted(bhn, key = lambda t:t[2])[-1]
		hp = bhn[1]
		ow.write(i + '\t' + hp + '\n')
	ow.close()

def compute_read_partitioning_accuracy(path, tagTXT):
	#path = os.path.abspath(sys.argv[1])
	subprocess.call('cat %s/bc1/*.allreads > %s/bc1/all.allreads'%(path, path), shell = True)
	allreads = '%s/bc1/all.allreads'%path
	#tagTXT = sys.argv[2]
	partitioned, tagged, readIsec = readBothFile(allreads, tagTXT)
	partitioned = partition_selection(partitioned, readIsec)
	print(len(partitioned), len(tagged), len(readIsec))
	pred_hap1, pred_hap2, considered_total, read_block_hp_nvar = readPred(partitioned, readIsec)
	trueHap1, trueHap2 = readTrue(tagged, readIsec)
	compute_accuracy(pred_hap1, pred_hap2, trueHap1, trueHap2, read_block_hp_nvar)

