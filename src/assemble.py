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

./assemble.py --illumina1 mom_1.fq dad_1.fq child_1.fq \\
              --illumina2 mom_2.fq dad_2.fq child_2.fq \\
              --pacbio pacbio_mom.fastq.gz

# Example2:

python assemble.py --illumina1 mom_1.fq dad_1.fq child_1.fq \\
                   --illumina2 mom_2.fq dad_2.fq child_2.fq \\
                   --pacbio pacbio_mom.fastq.gz pacbio_dad.fastq.gz pacbio_child.fastq.gz \\
                   --ped ped.ped

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
parser.add_argument('--illumina1', metavar='FASTQ/FASTA', type=str, nargs=3,
                     required = True, 
                     help='Illumina read data for right reads.')
parser.add_argument('--illumina2', metavar='FASTQ/FASTA', type=str, nargs=3,
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
parser.add_argument('-k', type = int, default = 77 , help = 'something')
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

logging.info('Concatenating fastq files for all left reads.')
for i in args.illumina1:
    logging.debug(i)
    c = subprocess.call(['cat %s >> %s/fq1.fq'%(i, tempPath)], shell = True)
logging.info('Concatenating fastq files for all right reads.')
for i in args.illumina2:
    logging.debug(i)
    c = subprocess.call(['cat %s >> %s/fq2.fq'%(i, tempPath)], shell = True)

# spades
logging.info('Running spades...')
spades_cmd = "python2 trioasm/SPAdes-3.13.0/spades.py -t 16 -k %d -1 %s/fq1.fq -2 %s/fq2.fq --only-assembler -o %s/illumina/" % (args.k, tempPath, tempPath, tempPath)
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
subprocess.call('%s snarls -l -t -r %s/illumina/asm1.trans %s/illumina/asm1.vg > %s/illumina/asm1.snarls' % (vg, tempPath, tempPath, tempPath), shell = True)


logging.info('Aligning...')
# alignment
if args.ped != None:
    for i in range(3):
        a = subprocess.call("%s -g %s/illumina/asm1.gfa -f %s -a %s/aln%d.gam --seeds-mum-count 100 --seeds-mxm-length 20" % (graphaligner, tempPath, args.pacbio[i], tempPath, i), shell = True, stdout=subprocess.PIPE)
        if a != 0:
            logging.error('Error while running GraphAligner. Exit Code:%d'%a)
            sys.exit(2)
elif args.ped == None:
    a = subprocess.call("%s -g %s/illumina/asm1.gfa -f %s -a %s/aln.gam --seeds-mum-count 100 --seeds-mxm-length 20" % (graphaligner, tempPath, args.pacbio[0], tempPath), shell = True, stdout=subprocess.PIPE)
    if a != 0:
        logging.error('Error while running GraphAligner. Exit Code: %d'%a)
        sys.exit(2) 

# global alignment
logging.info('Global alignment...')
# python3 issues not solved
if args.ped != None:
    for i in range(3):
        a = subprocess.call('python3 src/combine_canu_chunks_chunks.py %s/aln%d.gam %s/aln%d_g.gam'%(tempPath, i, tempPath, i), shell = True)
        if a != 0:
            logging.error('Error while running combine_canu_chunks_chunks.py. Exit Code:%d'%a)
            sys.exit(3)
elif args.ped == None:
    subprocess.call('python3 src/combine_canu_chunks_chunks.py %s/aln.gam %s/aln_g.gam'%(tempPath, tempPath), shell = True)
#subprocess.call('trioasm/wh_test/whatshap/venv/bin/python3 trioasm/combine_canu_chunks_chunks.py %s/alnc.gam %s/alnc_g.gam'%(tempPath, tempPath), shell = True)

logging.info('Bubble Chain...')
# bubble chain

if args.ped != None:
    subprocess.call("cat %s/aln0.gam %s/aln1.gam %s/aln2.gam > %s/aln.gam"%(tempPath, tempPath, tempPath, tempPath), shell = True)
subprocess.call('python3 src/find_bubble_chains_new_test.py %s/illumina/asm1.trans %s/aln.gam %s/bc1/aln.unitis > %s/bc1/aln.unitis.log 2>&1'%(tempPath, tempPath, tempPath, tempPath), shell = True)
subprocess.call('python3 src/aux_contigs_new_test.py %s/aln.gam %s/bc1/aln.unitis %s/bc1/aln.aux_unitis %s/illumina/asm1.trans'%(tempPath, tempPath, tempPath, tempPath), shell = True)
subprocess.call('python3 src/find_contigs_new.4.py %s/bc1/aln.unitis %s/bc1/aln.aux_unitis %s/bc1/aln.final_ctgs > %s/bc1/aln.final_ctgs.log'%(tempPath, tempPath, tempPath, tempPath), shell = True)
subprocess.call('python3 src/contigs_wh_test.py %s/illumina/asm1.trans %s/bc1/aln.final_ctgs 0 > %s/bc1/aln.log'%(tempPath, tempPath, tempPath), shell = True)

logging.info('Phaseg...')

# look for all .TRANS files
fileList = os.listdir(tempPath+ '/bc1')
translist = []
for i in fileList:
    if re.match(r'.*_contig.*.trans', i):
        translist.append(i)
print(translist)
if args.ped != None:
    for i in translist:
        #wh_test for 3 samples. Only 1 reads file
        
        subprocess.call("cd trioasm/whatshap_trioasm && python3 -m whatshap phaseg ../../%s/bc1/%s.reads ../../ped.ped ../../%s/bc1/%s ../../%s/aln0_g.gam ../../%s/aln1_g.gam ../../%s/aln2_g.gam > ../../%s.log" % (tempPath, i, tempPath, i, tempPath, tempPath, tempPath, i), shell = True)
elif args.ped == None:
    for i in translist:
        #MAV for single sample
        #subprocess.call('cd trioasm/whatshap', shell = True)
        subprocess.call('cd trioasm/whatshap && python3 -m whatshap phaseg ../../%s/bc1/%s.reads ../../%s/bc1/%s ../../%s/aln_g.gam ../../%s/illumina/asm1.vg > ../../%s.log'%(tempPath, i, tempPath, i, tempPath, tempPath, i), shell = True)#, args.reference, args.output), shell = True)
