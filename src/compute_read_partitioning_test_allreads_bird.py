# from collections import defaultdict
# import sys


# f = sys.argv[1]
# def compute_read_partitioning_accuracy():
# 	child_true_hap1 =[]
# 	child_true_hap2= []
# 	mom_true_hap1 =[]
# 	mom_true_hap2= []
# 	dad_true_hap1 =[]
# 	dad_true_hap2= []
# 	true_haps = defaultdict()
# 	#fromkeys(['child true_hap1', 'child true_hap2', 'mom true_hap1', 'dad true_hap1', 'dad true_hap2'])
# 	lines = []
# 	#for element in ['mom','dad','child']:
# 	for line in open(f):
# 		var = line.rstrip().split(" ")
# 		lines.append(var)
# 	for element in ['mom','dad','child']:
# 		for var in lines:
# 			if str(var[0][-1]) =='1':
# 				true_haps[element + 'true_hap1'].append(var[0])
# 			else:
# 				true_haps[element + 'true_hap2'].append(var[0])

# 	# pred_hap1 = defaultdict(list)
# 	# pred_hap2 = defaultdict(list)
# 	# total=0
# 	# blocks =set()
# 	# for line in open(f):
# 	# 	tokens= line.rstrip().split(" ")
# 	# 	blocks.add(tokens[1])
# 	# 	#print(tokens[1])
# 	# 	if tokens[4] =='1':
# 	# 		pred_hap1[tokens[4]].append(tokens[0])
# 	# 	else:
# 	# 		pred_hap2[tokens[4]].append(tokens[0])
# 	# 	total+=1

# 	# count = 0
# 	# #duplicates = [24899,18733,3763,5241,18841,17801,4623,17491,12705,1679,13699,18466,19418,24090,23744,13931,15366,13388,6334,2,8913,10192,20615,24999,24254,9050,15084,24061,5556,24708,10257,23275,18922,19628,23902,8623,22461,13405,21608,4279,17950,15274,15598]
# 	# for k in blocks:
# 	# #	if k in duplicates:
# 	# #		continue
# 	# 	count = count+max(len(list(set(pred_hap1[k]).intersection(set(true_hap1)))), len(list(set(pred_hap1[k]).intersection(set(true_hap2)))))
# 	# 	#print(k, count)
# 	# 	count= count+ max(len(list(set(pred_hap2[k]).intersection(set(true_hap1)))), len(list(set(pred_hap2[k]).intersection(set(true_hap2)))))
# 	# 	#print(k, count)
# 	# 	#print(k)
# 	# percent_partitionining_accuracy =  count/total
# 	# print('percent_partitionining_accuracy')
# 	# print(percent_partitionining_accuracy)

# compute_read_partitioning_accuracy()
from collections import defaultdict
import sys

f = sys.argv[1]
b3 = sys.argv[2]
b4 = sys.argv[3]

def compute_read_partitioning_accuracy():
	count = 0
	child_count = 0
	child_total = 0
	total=0


	true_hap1 = open(b3, 'r').read().split('\n')
	true_hap2 = open(b4, 'r').read().split('\n')
			
	pred_hap1 = defaultdict(list)
	pred_hap2 = defaultdict(list)
	blocks =set()
	for line in open(f):
		tokens = line.rstrip().split("\t")
		if int(tokens[1]) >= 50:
			try:
				total += 1
				blocks.add(tokens[3].split(':')[2])
				if tokens[2].split(':')[2] == '1':
					pred_hap1[tokens[1]].append(tokens[0])
				else:
					pred_hap2[tokens[1]].append(tokens[0])
			except IndexError:
				pass
	for k in blocks:
		#print('in count', count)
		#	if k in duplicates:
		#		continue
		count = count+max(len(list(set(pred_hap1[k]).intersection(set(true_hap1)))), len(list(set(pred_hap1[k]).intersection(set(true_hap2)))))
		#print(k, count)
		count = count+ max(len(list(set(pred_hap2[k]).intersection(set(true_hap1)))), len(list(set(pred_hap2[k]).intersection(set(true_hap2)))))


	percent_partitionining_accuracy =  count/total


compute_read_partitioning_accuracy()


