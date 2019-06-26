'''
Validate partitioning result for either simulated data or real data
'''
import sys, os
import re
import subprocess
from subprocess import PIPE
from . import validate_simulate
from . import validate_real

def valSim(path, pacbioFA):
    whdenovoPath = '/'.join(sys.path[0].split('/')[:-1])
    fileList = os.listdir(path)
    readslist = []
    allreadslist = []
    for i in fileList:
        if re.match(r'.*_block_*.allreads', i):
            allreadslist.append(i)
    subprocess.call('cat %s/bc1/*.allreads > %s/bc1/all.allreads'%(path, path), shell = True)
    validate_simulate.compute_read_partitioning_accuracy('%s/bc1/all.allreads'%path)
    #subprocess.call('python %s/src/validate_simulate.py %s/all.allreads'%(whdenovoPath, path), shell = True)
    #subprocess.call("cat %s/bc1/all.allreads | cut -d' ' -f1 | sort -u | uniq > %s/bc1/all.tmpreads"%(path, path), shell = True)
    #subprocess.call('python %s/src/get_unpartitioned.py %s/bc1/all.tmpreads %s'%(whdenovoPath, path, pacbioFA), shell = True)

def valReal(path, tag):
    validate_real.compute_read_partitioning_accuracy(path, tag)

def add_arguments(parser):
    arg = parser.add_argument
    arg('-p', '--path', required = True, help = 'The output path of partitioning results.')
    arg('-f', '--fasta', required = False, help = 'Required for simulated data validation.')
    arg('-t', '--truth', required = False, help = 'Required for real data validation.')

def main(args):
    if args.fasta != None and args.truth != None:
        print('choose only one!')
    elif args.fasta != None:
        valSim(args.path, args.fasta)
    else:
        valReal(args.path, args.truth)
    
