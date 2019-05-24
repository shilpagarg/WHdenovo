#! /usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
import subprocess
from subprocess import PIPE
import os, sys, time
import logging
import re

eg = '''
# Example1: 

./assemble.py --illumina1 mom_1.fq\\
              --illumina2 mom_2.fq\\
              --pacbio pacbio_mom.fastq.gz\\
              -s 100k -t 3

# Example2:

python assemble.py --illumina1 mom_1.fq dad_1.fq child_1.fq \\
                   --illumina2 mom_2.fq dad_2.fq child_2.fq \\
                   --pacbio pacbio_mom.fastq.gz pacbio_dad.fastq.gz pacbio_child.fastq.gz \\
                   --ped ped.ped -s 100k -t 3

# Order of MOM -> DAD -> CHILD matters!
'''

# prerequisite: stream, google, networkx, decorator, pysam

logging.basicConfig(level=logging.DEBUG)
logging.debug('Using python version %s'%sys.version)
if sys.version_info.major != 3:
    logging.error('Please use Python3.x to run this pipeline.')
    exit()
PID = os.getpid()
logging.info('PID %d'%PID)

parser = argparse.ArgumentParser(description=eg, formatter_class=RawTextHelpFormatter)
parser.add_argument('--illumina1', metavar='FASTQ/FASTA', type=str, nargs='+',
                     required = True, 
                     help='Illumina read data for right reads.')
parser.add_argument('--illumina2', metavar='FASTQ/FASTA', type=str, nargs='+',
                     required = True, 
                     help='Illumina read data for left reads.')
parser.add_argument('--pacbio', metavar='FASTQ/FASTA', type=str, nargs='+',
                     required = True, help='PacBio read data. If no PED file \n\
specified, should only give one file. \nOtherwise, three \
files in the order of \nMOM, DAD, CHILD.')
parser.add_argument('-p', '--ped', metavar='ped file', type=str,
                     required = False, help = "PED file.")
parser.add_argument('-o', '--output', metavar = 'PRED_HAPLOTIGS', 
                    required = False)
parser.add_argument('-f', '--reference', metavar = 'CANU', required = False)
#parser.add_argument('-s', '--size', required = True, help = 'Genome size.')
parser.add_argument('-k', type = int, default = 77 , help = 'something for Spades.')
parser.add_argument('-t', type = int, default = 3, help = "Some steps utilize GNU parallel, this is for specitying the threads for them.")
args = parser.parse_args()

if args.ped != None:
    if len(args.pacbio) != 3:
        raise ValueError('When PED specified, you should have three PACBIO '
                         'files in the order of MOM, DAD, CHILD.')
    if args.output != None or args.reference != None:
        raise ValueError('When PED is specified, you cannot have -o or -f.')
elif args.ped == None:
    if len(args.pacbio) != 1:
        raise TypeError('When no PED specified, you should have only one '
                        'PACBIO file.')
    #if args.output == None or args.reference == None:
    #    raise ValueError('When we don\'t have PED, we should have both -o and -f')

# TODO check if file amount match between illumina1 illumina2 pacbio. 


vg = "trioasm/vg"
graphaligner = "trioasm/Aligner"
# spades commands
tempPath = 'temp_%d_%s'%(PID, time.strftime("%m%d%H%M%S"))

logging.info('Making temp directories.')
if not os.path.isdir(tempPath):
    subprocess.call(['mkdir', tempPath])
    subprocess.call(['mkdir', '%s/illumina'%tempPath])
    subprocess.call(['mkdir', '%s/bc1'%tempPath])
if os.path.exists('%s/fq1.fq'%tempPath):
    subprocess.call(['rm', '%s/fq1.fq'%tempPath])
if os.path.exists('%s/fq2.fq'%tempPath):
    subprocess.call(['rm', '%s/fq2.fq'%tempPath])


for i in args.illumina1:
    logging.debug(i)
    c = subprocess.call(['cat %s >> %s/fq1.fq'%(i, tempPath)], shell = True)
logging.info('Concatenating fastq files for all right reads.')
for i in args.illumina2:
    logging.debug(i)
    c = subprocess.call(['cat %s >> %s/fq2.fq'%(i, tempPath)], shell = True)

# spades  ADDED "-m 500" FOR LARGE GENOME
logging.info('Running spades...')
spades_cmd = "python2 trioasm/SPAdes-3.13.0/spades.py -t %d -k %d -m 500 -1 %s/fq1.fq -2 %s/fq2.fq --only-assembler -o %s/illumina/" % (args.t, args.k, tempPath, tempPath, tempPath)
spades_cmd = spades_cmd.split()
a = subprocess.call(spades_cmd, shell=False)#, stdout=subprocess.PIPE)
if a != 0:
    logging.error('Error while running spades. Error Code: %d'%a)
    sys.exit(1)

# filter graph
logging.info('Filtering graph...')
subprocess.call("grep -v '^P' %s/illumina/assembly_graph_with_scaffolds.gfa | awk -F'\\t' '{ if ($2 != $4) print $0}' | %s view --gfa-in - --vg | %s view -g - | awk -F'\\t' '{ if ($2 !=$4) print $0}' > %s/asm1.gfa" % (tempPath , vg, vg, tempPath), shell = True)
subprocess.call("python2 src/printnodedegrees_gfa.py %s/asm1.gfa | awk -F' ' '{ if($2 > 70 || $2==0) printf \"%%s\\n\", $1 }' > %s/asm1.wrongnodes"%(tempPath, tempPath), shell = True)
subprocess.call('python2 src/remove_wrongnodes.py %s/asm1.wrongnodes %s/asm1.gfa %s/illumina/asm1.gfa'%(tempPath, tempPath, tempPath), shell = True)

logging.info('Running snarls...')
# snarl
subprocess.call('%s view --gfa-in --vg %s/illumina/asm1.gfa > %s/illumina/asm1.vg' % (vg, tempPath, tempPath), shell = True)
subprocess.call('%s snarls -t -r %s/illumina/asm1.trans %s/illumina/asm1.vg > %s/illumina/asm1.snarls' % (vg, tempPath, tempPath, tempPath), shell = True)


logging.info('Aligning...')
# alignment
if args.ped != None:
    for i in range(3):
        a = subprocess.call("%s -t %d -g %s/illumina/asm1.gfa -f %s -a %s/aln%d.gam --seeds-mum-count 200 --seeds-mxm-length 20" % (graphaligner, args.t, tempPath, args.pacbio[i], tempPath, i), shell = True, stdout=subprocess.PIPE)
        if a != 0:
            logging.error('Error while running GraphAligner. Exit Code:%d'%a)
            sys.exit(2)
elif args.ped == None:
    a = subprocess.call("%s -t %d -g %s/illumina/asm1.gfa -f %s -a %s/aln.gam --seeds-mum-count 200 --seeds-mxm-length 20" % (graphaligner, args.t, tempPath, args.pacbio[0], tempPath), shell = True, stdout=subprocess.PIPE)
    if a != 0:
        logging.error('Error while running GraphAligner. Exit Code: %d'%a)
        sys.exit(2) 


logging.info('Bubble Chain...')

# bubble chain
# ========================================For 4scripts version, change the script name (bubble_chian.py) to bubble_chain_4scripts_version.py====================================
if args.ped != None:
    subprocess.call('python3 src/bubble_chain_directed.py %s/illumina/asm1.trans %s/aln0.gam %s/aln1.gam %s/aln2.gam %s/bc1/aln'%(tempPath, tempPath, tempPath, tempPath, tempPath), shell = True)
    #subprocess.call('python3 src/bubble_chain_4scripts_version.py %s/illumina/asm1.trans %s/aln0.gam %s/aln1.gam %s/aln2.gam %s/bc1/aln'%(tempPath, tempPath, tempPath, tempPath, tempPath), shell = True)
    for i in range(3):
        subprocess.call("ls %s/bc1/*trans | sed 's/.trans//' | parallel -j %d 'python3 src/combine_canu_chunks_chunks.py {}_%d.gam {}_%d_g.gam'" % (tempPath, args.t, i, i), shell = True)
else:
    subprocess.call('python3 src/bubble_chain_directed.py %s/illumina/asm1.trans %s/aln.gam %s/bc1/aln'%(tempPath, tempPath, tempPath), shell = True)
    subprocess.call("ls %s/bc1/*trans | sed 's/.trans//' | parallel -j %d 'python3 src/combine_canu_chunks_chunks.py {}.gam {}_g.gam'" % (tempPath, args.t), shell = True)


logging.info('Phaseg...')

if args.ped != None:
    subprocess.call("cd trioasm/whatshap_trioasm && ls ../../%s/bc1/*trans | sed 's/.trans//' | parallel -j %d 'python3 -m whatshap phaseg {}.reads ../../ped.ped {}.trans {}_0_g.gam {}_1_g.gam {}_2_g.gam'" % (tempPath, args.t), shell = True)
elif args.ped == None:
    subprocess.call("cd trioasm/whatshap && ls ../../%s/bc1/*trans | sed 's/.trans//' | parallel -j %d 'python3 -m whatshap phaseg {}.reads {}.trans {}_g.gam ../../%s/illumina/asm1.vg > {}.log'" % (tempPath, args.t, tempPath), shell = True)#, args.reference, args.output), shell = True)
