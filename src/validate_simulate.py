from collections import defaultdict
import sys

f = sys.argv[1]
def compute_read_partitioning_accuracy():
	count = 0
	child_count = 0
	child_total = 0
	total=0
	for element in ['mom','dad','child']:
		true_hap1 =[]
		true_hap2= []
		
		for line in open(f):
			var = line.rstrip().split(" ")
			if element in var[0]:
				if 'child1' in str(var[0]) or 'dad1' in str(var[0]) or 'mom1' in str(var[0]):
					true_hap1.append(var[0])
				else:
					true_hap2.append(var[0])
			
		pred_hap1 = defaultdict(list)
		pred_hap2 = defaultdict(list)
		blocks =set()
		for line in open(f):
			tokens = line.rstrip().split(" ")
			if element in tokens[0]:
				blocks.add(tokens[1])
				if tokens[2] =='1':
					pred_hap1[tokens[1]].append(tokens[0])
				else:
					pred_hap2[tokens[1]].append(tokens[0])
				total += 1
				if element == 'child':
					child_total += 1
		for k in blocks:
			count += max(len(list(set(pred_hap1[k]).intersection(set(true_hap1)))), len(list(set(pred_hap1[k]).intersection(set(true_hap2)))))
			count += max(len(list(set(pred_hap2[k]).intersection(set(true_hap1)))), len(list(set(pred_hap2[k]).intersection(set(true_hap2)))))

			if element == 'child':
				child_count += max(len(list(set(pred_hap1[k]).intersection(set(true_hap1)))), len(list(set(pred_hap1[k]).intersection(set(true_hap2)))))
				child_count += max(len(list(set(pred_hap2[k]).intersection(set(true_hap1)))), len(list(set(pred_hap2[k]).intersection(set(true_hap2)))))
	print('final count',count)
	print('total', total)
	print('child_count', child_count)
	print('child_total', child_total)
	child_partitioning_accuracy = child_count/child_total
	print('child_partitionining_accuracy')
	print(child_partitioning_accuracy)
	percent_partitionining_accuracy =  count/total
	print('percent_partitionining_accuracy')
	print(percent_partitionining_accuracy)

compute_read_partitioning_accuracy()


