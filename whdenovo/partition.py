'''
From input short- and long- reads, run the whole partitioning pipeline and get the final partitioned read list.
'''
import argparse
from argparse import RawTextHelpFormatter
import subprocess
from subprocess import PIPE
import os, sys, time
import logging
import re

def checkStatus(args):
    if sys.version_info.major != 3:
        logging.error('Please use Python3.x to run this pipeline.')
        exit()
    if args.ped != None:
        if len(args.pacbio) != 3:
            raise ValueError('When PED specified, you should have three PACBIO '
                             'files in the order of MOM, DAD, CHILD.')
    elif args.ped == None:
        if len(args.pacbio) != 1:
            raise ValueError('When no PED specified, you should have only one '
                            'PACBIO file.')
    if not os.path.exists(args.illumina1):
        raise FileNotFoundError('No such file or directory: %s' % args.illumina1)
    if not os.path.exists(args.illumina2):
        raise FileNotFoundError('No such file or directory: %s' % args.illumina2)
    for file in args.pacbio:
        if not os.path.exists(file):
            raise FileNotFoundError('No such file or directory: %s'%file)

def run_pipeline(args):

    logging.basicConfig(level=logging.DEBUG)
    logging.info('Using python version %s'%sys.version)
    
    PID = os.getpid()
    logging.info('PID %d'%PID)
    
    whdenovoPath = '/'.join(sys.path[0].split('/')[:-1])
    whdenovoPath = sys.path[0]
    vg = "%s/trioasm/vg"%whdenovoPath
    logging.info('Using vg at: %s'%vg)
    graphaligner = "GraphAligner"
    if args.output != None:
        tempPath = os.getcwd() + '/' + args.output
    else:
        tempPath = os.getcwd() + '/temp_%d_%s'%(PID, time.strftime("%m%d%H%M%S"))

    logging.info('Making temp directories.')
    if not os.path.isdir(tempPath):
        subprocess.call(['mkdir', tempPath])
        subprocess.call(['mkdir', '%s/illumina'%tempPath])
        subprocess.call(['mkdir', '%s/bc1'%tempPath])

    # spades  ADDED "-m 500" FOR LARGE GENOME
    #logging.info('bfc error correcting...')
    #subprocess.call('bfc -t %d -s %s %s 1> %s/cor.1.fq 2>> %s/bfc.log'%(args.t, args.size, args.illumina1, tempPath, tempPath), shell = True)
    #subprocess.call('bfc -t %d -s %s %s 1> %s/cor.2.fq 2>> %s/bfc.log'%(args.t, args.size, args.illumina2, tempPath, tempPath), shell = True)
    #logging.info('bfc log saved at %s/bfc.log'%tempPath)
     
    logging.info('Running spades...')
    spades_cmd = "python2 %s/trioasm/SPAdes-3.13.0/spades.py -t %d -k %d -m 500 -1 %s/cor.1.fq -2 %s/cor.2.fq --only-assembler -o %s/illumina/" % (whdenovoPath, args.t, args.k, tempPath, tempPath, tempPath)
    spades_cmd = spades_cmd.split()
    a = subprocess.call(spades_cmd, shell=False, stdout=subprocess.PIPE)
    if a != 0:
        logging.error('Error while running spades. Error Code: %d'%a)
        sys.exit(1)
    logging.info('SPAdes log saved at %s/illumina/spades.log'%tempPath)

    logging.info('Filtering graph...')
    subprocess.call("grep -v '^P' %s/illumina/assembly_graph_with_scaffolds.gfa | awk -F'\\t' '{ if ($2 != $4) print $0}' | %s view --gfa-in - --vg | %s view -g - | awk -F'\\t' '{ if ($2 !=$4) print $0}' > %s/asm1.gfa" % (tempPath , vg, vg, tempPath), shell = True)
    subprocess.call("python2 %s/whdenovo/printnodedegrees_gfa.py %s/asm1.gfa | awk -F' ' '{ if($2 > 70 || $2==0) printf \"%%s\\n\", $1 }' > %s/asm1.wrongnodes"%(whdenovoPath, tempPath, tempPath), shell = True)
    subprocess.call('python2 %s/whdenovo/remove_wrongnodes.py %s/asm1.wrongnodes %s/asm1.gfa %s/illumina/asm1.gfa'%(whdenovoPath, tempPath, tempPath, tempPath), shell = True)
    logging.info('Running snarls...')
    subprocess.call('%s view --gfa-in --vg %s/illumina/asm1.gfa > %s/illumina/asm1.vg' % (vg, tempPath, tempPath), shell = True)
    subprocess.call('%s snarls -t -r %s/illumina/asm1.trans %s/illumina/asm1.vg > %s/illumina/asm1.snarls' % (vg, tempPath, tempPath, tempPath), shell = True)

    logging.info('Aligning...')
    if args.ped != None:
        for i in range(3):
            a = subprocess.call("%s -t %d -g %s/illumina/asm1.gfa -f %s -a %s/aln%d.gam --seeds-mum-count 100000 --seeds-mxm-length 10 -C 500000 -b 35" % (graphaligner, args.t, tempPath, args.pacbio[i], tempPath, i), shell = True, stdout=subprocess.PIPE)
            if a != 0:
                logging.error('Error while running GraphAligner. Exit Code:%d'%a)
                sys.exit(2)
    elif args.ped == None:
        a = subprocess.call("%s -t %d -g %s/illumina/asm1.gfa -f %s -a %s/aln.gam --seeds-mum-count 100000 --seeds-mxm-length 10 -C 500000 -b 35" % (graphaligner, args.t, tempPath, args.pacbio[0], tempPath), shell = True, stdout=subprocess.PIPE)
        if a != 0:
            logging.error('Error while running GraphAligner. Exit Code: %d'%a)
            sys.exit(a) 
    a = subprocess.call("ls %s/aln*gam | parallel '%s view -a {} > {}.json'"%(tempPath, vg), shell = True)
    if a != 0:
        logging.error('Error while converting GAM to JSON. Exit Code: %d'%a)
        sys.exit(a)

    logging.info('Partitioning...')
    if args.ped != None:
        subprocess.call("cd %s/trioasm/whatshap_trioasm && python -m whatshap phaseg reads %s %s/illumina/asm1.trans %s/aln0.gam.json %s/aln1.gam.json %s/aln2.gam.json -t %d -p %s/bc1/aln --lowc %d --high %d > %s/partition.log" % (whdenovoPath, args.ped, tempPath, tempPath, tempPath, tempPath, args.t, tempPath, args.lowc, args.highc, tempPath), shell = True)
        subprocess.call("ls %s/bc1/*allreads | parallel -j%d \"awk '\\$3 == 1 {print \\$1}' {} > {}.hp1.reads\"" % (tempPath, args.t), shell = True)
        subprocess.call("ls %s/bc1/*allreads | parallel -j%d \"awk '\\$3 == 0 {print \\$1}' {} > {}.hp0.reads\"" % (tempPath, args.t), shell = True)
        subprocess.call("cat %s/bc1/*hp1.reads | sort | uniq > %s/HP1.reads" % (tempPath,tempPath), shell = True)
        subprocess.call("cat %s/bc1/*hp0.reads | sort | uniq > %s/HP0.reads" % (tempPath,tempPath), shell = True)
        logging.info('Partitioning finished. Read names of different haplotypes are saved in:')
        print('%s/HP0.reads'%tempPath)
        print('%s/HP1.reads'%tempPath)
    else:
        # TODO individual case
        pass
        
def run_test(sth):
    print(sth)

def add_arguments(parser):
    arg = parser.add_argument
    arg('--illumina1', metavar='FASTQ/FASTA', type=str, 
                         required = True, 
                         help='Illumina short-read data for right reads.')
    arg('--illumina2', metavar='FASTQ/FASTA', type=str, 
                         required = True, 
                         help='Illumina short-read data for left reads.')
    arg('--pacbio', metavar='FASTQ/FASTA', type=str, nargs='+',
                         required = True, help='PacBio long-read data. If no PED file specified, should only give one file. \nOtherwise, three \
    files in the order of MOM, DAD, CHILD.')
    arg('-p', '--ped', metavar='ped file', type=str,
                         required = False, help = "PED file.")
    arg('-o', '--output', metavar = 'PATH', 
                        required = False)
    arg('-k', type = int, default = 77 , help = 'K-mer size setting for SPAdes, must be odd and no larger than 128. [77].')
    arg('-s', '--size', type = str, required = True, help = 'Expected genome size, acceptible example: 50k, 24m, 2g.')
    arg('-t', metavar = 'INT', type = int, default = 4, help = "Use multiprocessing in the algorithm, and some steps utilize GNU parallel. [4]")
    arg('--lowc', metavar = 'INT', type = int, default = 5, help = 'Lowest threshold for coverage to support edges.')
    arg('--highc', metavar = 'INT', type = int, default = 20, help = 'Highest threshold for coverage to detect repeats.')

def main(args):
    checkStatus(args)
    run_pipeline(args)
    
