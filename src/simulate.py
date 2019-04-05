#! /usr/bin/env python3
import sys, os
import subprocess
from subprocess import PIPE
import argparse

'''
IMPORTANT: Be sure your input original fasta file is in UNIX fileformat, 
           Otherwise the PBSIM program cannot interpret things properly.
'''
parser = argparse.ArgumentParser(description='Use this tool to simulate read data from a given FASTA sequence.')
sub = parser.add_subparsers(dest = 'subparser')

illumina = sub.add_parser('illumina')
illumina.add_argument('fasta', type=str, 
                     help='Original FASTA sequence file.')
illumina.add_argument('het', type=str, help = 'A NUMBER of het level.')
illumina.add_argument('dir', type = str, 
					  help = 'Output directory name.')

pacbio = sub.add_parser('pacbio', help = 'Take simulated ')
pacbio.add_argument('fastq', type = str, 
					help = 'A FASTQ read file.')
pacbio.add_argument('coverage', type = int, 
					help = 'A NUMBER for the coverage in the pacbio read output. ')
pacbio.add_argument('dir', type = str, 
					help = 'Output directory name.')
pacbio.add_argument('fasta', type = str, nargs = 6,
					help = 'Original FASTA sequence files. Please input them in the order of: MOM1 MOM2 DAD1 DAD2 CHILD1 CHILD2')

args = parser.parse_args()

family = ['mom1', 'mom2', 'dad1', 'dad2', 'child1', 'child2']
# ILLUMINA part
if args.subparser == 'illumina':

	mom1FA = args.fasta
	het = args.het
	illuminaOutput = os.path.abspath(args.dir)
	if not os.path.exists(mom1FA):
		raise FileNotFoundError(mom1FA + ' not found.')
		sys.exit(1)
	# Generate haplotypes.
	simulate_mut = 'src/simulate_mutation.py'
	subprocess.run('mkdir %s'%illuminaOutput, shell = True	)
	subprocess.call('python3 %s %s %s %s/mom2.het%s.cov30'%(simulate_mut, mom1FA, het, illuminaOutput, het), shell = True)
	subprocess.call('python3 %s %s %s %s/dad1.het%s.cov30'%(simulate_mut, mom1FA, het, illuminaOutput, het), shell = True)
	subprocess.call('python3 %s %s %s %s/dad2.het%s.cov30'%(simulate_mut, mom1FA, het, illuminaOutput, het), shell = True)
	subprocess.call('cp %s/mom2.het%s.cov30.fasta %s/child1.het%s.cov30.fasta'%(illuminaOutput, het, illuminaOutput, het), shell = True)
	subprocess.call('cp %s/dad1.het%s.cov30.fasta %s/child2.het%s.cov30.fasta'%(illuminaOutput, het, illuminaOutput, het), shell = True)
	subprocess.call('cp %s %s/mom1.het%s.cov30.fasta'%(mom1FA, illuminaOutput, het), shell = True)

	# simulate illumina reads
	art = 'trioasm/art_bin_MountRainier/art_illumina'
	for i in family:
		subprocess.call('%s  -ss HSXn -sam -i %s/%s.het%s.cov30.fasta -rs 777 -p -l 150 -f 15 -m 200 -s 10 -o %s/%s.het%s.cov30_'%(art, illuminaOutput, i, het, illuminaOutput, i, het), shell = True)

	# Combine left and right reads per individual.
	for i in ['mom', 'dad', 'child']:
		subprocess.call('cat %s/%s1.het%s.cov30_1.fq %s/%s2.het%s.cov30_1.fq > %s/%s.het%s.cov30_1.fq'%(illuminaOutput, i, het, illuminaOutput, i, het, illuminaOutput, i, het), shell = True)
		subprocess.call('cat %s/%s1.het%s.cov30_2.fq %s/%s2.het%s.cov30_2.fq > %s/%s.het%s.cov30_2.fq'%(illuminaOutput, i, het, illuminaOutput, i, het, illuminaOutput, i, het), shell = True)

# PACBIO part
elif args.subparser == 'pacbio':
	FAs = args.fasta
	fq = args.fastq
	halfCov = args.coverage / 2
	outputDir = os.path.abspath(args.dir)
	subprocess.run('mkdir ' + outputDir, shell = True)
	for i in range(6):
		a = subprocess.call('pbsim --seed 10 --prefix %s/%s --depth %d --sample-fastq %s %s > pbsim.log 2>&1'%(outputDir, family[i], halfCov, fq, FAs[i]), shell = True)
		subprocess.call("awk 'NR%%4==1 {printf(\"%%s_%s\\n\",$0)} NR%%4!=1 {print}\\' %s/%s_0001.fastq > %s/pacbio_%s.fastq" % (family[i], outputDir, family[i], outputDir, family[i]), shell = True)
		print('%s/pacbio_%s.fastq'%(outputDir, family[i]))
		if a != 0:
			print('Check if you have installed package package "pbsim". Then check if anything wrong with your input file.')
			sys.exit(2)
	for i in ['mom', 'dad', 'child']:
		subprocess.call('cat %s/pacbio_%s1.fastq %s/pacbio_%s2.fastq | seqtk seq -A - > %s/pacbio_%s.fasta'%(outputDir, i, outputDir, i, outputDir, i), shell = True)
		print('%s/pacbio_%s.fasta'%(outputDir, i))
	print('Above are all you need.')