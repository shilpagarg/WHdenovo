#! /usr/bin/env python3
import argparse
import subprocess

parser = argparse.ArgumentParser(description='whdenovo. You can simulate data, do assembly, and validate your assembly.')
sub = parser.add_subparsers(dest = 'sub')

assemble = sub.add_parser('assemble')
assemble.add_argument('--illumina1', metavar='FASTQ/FASTA', type=str, nargs=3,
                     required = True, 
                     help='Illumina read data for right reads.')
assemble.add_argument('--illumina2', metavar='FASTQ/FASTA', type=str, nargs=3,
                     required = True, 
                     help='Illumina read data for left reads.')
assemble.add_argument('--pacbio', metavar='FASTQ/FASTA', type=str, nargs='+',
                     required = True, help='PacBio read data. If no PED file \n\
specified, should only give one file. \nOtherwise, three \
files in the order of \nMOM, DAD, CHILD.')
assemble.add_argument('-p', '--ped', metavar='ped file', type=str,
                     required = False, help = "PED file.")
assemble.add_argument('-o', '--output', metavar = 'PRED_HAPLOTIGS', 
                    required = False)
assemble.add_argument('-f', '--reference', metavar = 'CANU', required = False)
assemble.add_argument('-k', type = int, default = 77 , help = 'something')

validate = sub.add_parser('validate')
validate.add_argument('outputDir', type = str, 
	                  help = 'e.g. temp_PID_TIME/bc1')
validate.add_argument('pacbio', type = str, 
	                  help = 'Child\'s PacBio FASTA reads.')

simulate = sub.add_parser('simulate')
simusub = simulate.add_subparsers(dest = 'simusub')

illumina = simusub.add_parser('illumina', help = 'Simulate Illumina reads of MOM1, MOM2, DAD1, DAD2, CHILD1, CHILD2.')
illumina.add_argument('fasta', type=str, 
                     help='Original FASTA sequence file.')
illumina.add_argument('het', type=str, help = 'A NUMBER of het level.')
illumina.add_argument('dir', type = str, 
					  help = 'Output directory name.')

pacbio = simusub.add_parser('pacbio', help = 'Take simulated haplotypes from "illumina" subcommand to simulate PacBio reads.')
pacbio.add_argument('fastq', type = str, 
					help = 'A FASTQ read file.')
pacbio.add_argument('coverage', type = int, 
					help = 'A NUMBER for the coverage in the pacbio read output. ')
pacbio.add_argument('dir', type = str, 
					help = 'Output directory name.')
pacbio.add_argument('fasta', type = str, nargs = 6,
					help = 'Original FASTA sequence files. Please input them in the order of: MOM1 MOM2 DAD1 DAD2 CHILD1 CHILD2')

args = parser.parse_args()

if args.sub == 'simulate':
	if args.simusub == 'illumina':
		subprocess.call('src/simulator.py illumina %s %s %s'%(args.fasta, args.het, args. dir), shell = True)
	elif args.simusub == 'pacbio':
		subprocess.call('src/simulator.py pacbio %s %s %s %s'%(args.fastq, args.coverage, args.dir, ' '.join(args.fasta)), shell = True)#[0], args.fasta[1], args.fasta[2], args.fasta[3], args.fasta[4], args.fasta[5]), shell = True)
elif args.sub == 'validate':
	subprocess.call('./validator.py %s %s'%(args.outputDir, args.pacbio), shell = True)
elif args.sub == 'assembly':
	if args.ped == None:
		subprocess.call('./assembly.py --illumina1 %s --illumina2 %s --pacbio %s'%(' '.join(args.illumina1), ' '.join(args.illumina2), args.pacbio[0]), shell = True)
	elif args.ped != None:
		subprocess.call('./assembly.py --illumina1 %s --illumina2 %s --pacbio %s --ped %s'%(' '.join(args.illumina1), ' '.join(args.illumina2), ' '.join(args.pacbio[0]), args.ped), shell = True)
