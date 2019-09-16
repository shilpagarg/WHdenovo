'''
Simulate Illumina short-read data, PacBio long-read data, or PacBioCCS long-accurate-read data.
'''
import sys, os
import subprocess
from subprocess import PIPE
import argparse
from . import simulate_mutation
'''
IMPORTANT: Be sure your input original fasta file is in UNIX fileformat, 
           Otherwise the PBSIM program cannot interpret things properly.
'''
helpMessage = '''usage: whdenovo simulate [-h] {illumina,pacbio,pacbioccs} ...

Simulate Illumina short-read data, PacBio long-read data, or PacBioCCS long-accurate-read data.

positional arguments:
  {illumina,pacbio,pacbioccs}
    illumina            Simulate haplotype sequences for three individuals
    pacbio              Take simulated illumina data and generate pacbio reads
    pacbioccs           Simulate haplotypes and corresponding pacbioccs reads

optional arguments:
  -h, --help            show this help message and exit'''
def sim(sub, args):
    
    pwd = os.getcwd()
    whdenovoPath = '/'.join(os.path.abspath(__file__).split('/')[:-1])
    family = ['mom1', 'mom2', 'dad1', 'dad2', 'child1', 'child2']
    # ILLUMINA part
    if sub == 'illumina':

        mom1FA = args.fasta
        het = args.het
        illuminaOutput = pwd + '/' + args.dir
        subprocess.run('mkdir -p ' + illuminaOutput, shell = True)
        if not os.path.exists(mom1FA):
            raise FileNotFoundError(mom1FA + ' not found.')
            sys.exit(1)
    	# Generate haplotypes.
        simulate_mutation.simMut(mom1FA, het, '%s/mom2.het%s.cov30' % (illuminaOutput, het))
        simulate_mutation.simMut(mom1FA, het, '%s/dad1.het%s.cov30' % (illuminaOutput, het))
        simulate_mutation.simMut(mom1FA, het, '%s/dad2.het%s.cov30' % (illuminaOutput, het))
        #simulate_mut = whdenovoPath + '/simulate_mutation.py'
        #subprocess.call('python3 %s %s %s %s/mom2.het%s.cov30'%(simulate_mut, mom1FA, het, illuminaOutput, het), shell = True)
        #subprocess.call('python3 %s %s %s %s/dad1.het%s.cov30'%(simulate_mut, mom1FA, het, illuminaOutput, het), shell = True)
        #subprocess.call('python3 %s %s %s %s/dad2.het%s.cov30'%(simulate_mut, mom1FA, het, illuminaOutput, het), shell = True)
        subprocess.call('cp %s/mom2.het%s.cov30.fasta %s/child1.het%s.cov30.fasta'%(illuminaOutput, het, illuminaOutput, het), shell = True)
        subprocess.call('cp %s/dad1.het%s.cov30.fasta %s/child2.het%s.cov30.fasta'%(illuminaOutput, het, illuminaOutput, het), shell = True)
        subprocess.call('cp %s %s/mom1.het%s.cov30.fasta'%(mom1FA, illuminaOutput, het), shell = True)

        # simulate illumina reads
        art = 'art_illumina'
        for i in family:
            subprocess.call('%s  -ss HSXn -sam -i %s/%s.het%s.cov30.fasta -p -l 150 -f 15 -m 200 -s 10 -o %s/%s.het%s.cov30_'%(art, illuminaOutput, i, het, illuminaOutput, i, het), shell = True, stdout=subprocess.PIPE)

    	# Combine left and right reads per individual.
        for i in ['mom', 'dad', 'child']:
            subprocess.call('cat %s/%s1.het%s.cov30_1.fq %s/%s2.het%s.cov30_1.fq > %s/%s.het%s.cov30_1.fq'%(illuminaOutput, i, het, illuminaOutput, i, het, illuminaOutput, i, het), shell = True)
            subprocess.call('cat %s/%s1.het%s.cov30_2.fq %s/%s2.het%s.cov30_2.fq > %s/%s.het%s.cov30_2.fq'%(illuminaOutput, i, het, illuminaOutput, i, het, illuminaOutput, i, het), shell = True)
        print('WHAT YOU NEED FOR PACBIO SIMULATION: ')
        print('%s/mom1.het%s.cov30.fasta %s/mom2.het%s.cov30.fasta %s/dad1.het%s.cov30.fasta %s/dad2.het%s.cov30.fasta %s/child1.het%s.cov30.fasta %s/child2.het%s.cov30.fasta' % (illuminaOutput, het, illuminaOutput, het, illuminaOutput, het, illuminaOutput, het, illuminaOutput, het, illuminaOutput, het))
        print('WHAT YOU NEED FOR ASSEMBLING: ')
        print('%s/mom.het%s.cov30_1.fq\t%s/dad.het%s.cov30_1.fq\t%s/child.het%s.cov30_1.fq'%(illuminaOutput, het, illuminaOutput, het, illuminaOutput, het))
        print('%s/mom.het%s.cov30_2.fq\t%s/dad.het%s.cov30_2.fq\t%s/child.het%s.cov30_2.fq'%(illuminaOutput, het, illuminaOutput, het, illuminaOutput, het))

    # PACBIO part
    elif sub == 'pacbio':
        FAs = args.fasta
        fq = args.fastq
        halfCov = args.coverage / 2
        outputDir = pwd + '/' + args.dir
        subprocess.run('mkdir -p ' + outputDir, shell = True)
        for i in range(6):
            a = subprocess.call('pbsim --seed %d --prefix %s/%s --depth %d --sample-fastq %s %s > %s/pbsim.log 2>&1'%(i+1, outputDir, family[i], halfCov, fq, FAs[i], outputDir), shell = True)
            subprocess.call("awk 'NR%%4==1 {printf(\"%%s_%s\\n\",$0)} NR%%4!=1 {print}\\' %s/%s_0001.fastq > %s/pacbio_%s.fastq" % (family[i], outputDir, family[i], outputDir, family[i]), shell = True)
            print('%s/pacbio_%s.fastq'%(outputDir, family[i]))
            if a != 0:
                print('Check if you have installed package package "pbsim". Then check if anything wrong with your input file.')
                sys.exit(2)
        for i in ['mom', 'dad', 'child']:
            subprocess.call('cat %s/pacbio_%s1.fastq %s/pacbio_%s2.fastq | seqtk seq -A - > %s/pacbio_%s.fasta'%(outputDir, i, outputDir, i, outputDir, i), shell = True)
            print('%s/pacbio_%s.fasta'%(outputDir, i))
        print('Above are all you need.')

    elif sub == 'pacbioccs':
            mom1FA = args.fasta
            fq = 'test/pacbioccs.sample.fastq'
            het = args.het
            pbccsOutput = pwd + '/' + args.dir
            subprocess.run('mkdir -p ' + pbccsOutput, shell = True)
            if not os.path.exists(mom1FA):
                    raise FileNotFoundError(mom1FA + ' not found.')
                    sys.exit(1)
            # Generate haplotypes.
            simulate_mutation.simMut(mom1FA, het, '%s/mom2.het%s.cov30' % (pbccsOutput, het))
            simulate_mutation.simMut(mom1FA, het, '%s/dad1.het%s.cov30' % (pbccsOutput, het))
            simulate_mutation.simMut(mom1FA, het, '%s/dad2.het%s.cov30' % (pbccsOutput, het))    
            #simulate_mut = whdenovoPath + '/simulate_mutation.py'
            #subprocess.call('python3 %s %s %s %s/mom2'%(simulate_mut, mom1FA, het, pbccsOutput), shell = True)
            #subprocess.call('python3 %s %s %s %s/dad1'%(simulate_mut, mom1FA, het, pbccsOutput), shell = True)
            #subprocess.call('python3 %s %s %s %s/dad2'%(simulate_mut, mom1FA, het, pbccsOutput), shell = True)
            subprocess.call('cp %s/mom2.fasta %s/child1.fasta'%(pbccsOutput, pbccsOutput), shell = True)
            subprocess.call('cp %s/dad1.fasta %s/child2.fasta'%(pbccsOutput, pbccsOutput), shell = True)
            subprocess.call('cp %s %s/mom1.fasta'%(mom1FA, pbccsOutput), shell = True)
            FAs = ['mom1.fasta',  'mom2.fasta', 'dad1.fasta', 'dad2.fasta', 'child1.fasta', 'child2.fasta']

            for i in range(6):
                    a = subprocess.call('pbsim --seed %d --prefix %s/%s --depth 10 --sample-fastq %s %s/%s > %s/pbsim.log 2>&1'%(i, pbccsOutput, family[i], fq, pbccsOutput, FAs[i], pbccsOutput), shell = True
    )
                    subprocess.call("awk 'NR%%4==1 {printf(\"%%s_%s\\n\",$0)} NR%%4!=1 {print}\\' %s/%s_0001.fastq > %s/pacbioccs_%s.fastq" % (family[i], pbccsOutput, family[i], pbccsOutput, family[i]
    ), shell = True)
                    #print('%s/pacbioccs_%s.fastq'%(outputDir, family[i]))
                    if a != 0:
                            print('Check if you have installed package package "pbsim". Then check if anything wrong with your input file.')
                            sys.exit(2)
            for i in ['mom', 'dad', 'child']:
                    subprocess.call('cat %s/pacbioccs_%s1.fastq %s/pacbioccs_%s2.fastq | seqtk seq -A - > %s/pacbioccs_%s.fasta'%(pbccsOutput, i, pbccsOutput, i, pbccsOutput, i), shell = True)
                    print('%s/pacbioccs_%s.fasta'%(pbccsOutput, i))
            print('Above are all you need.')

def add_arguments(parser):
    #arg = parser.add_argument
    #arg('-p', required = True, help = 'print this string')

    sub = parser.add_subparsers(dest = 'subparser')

    illumina = sub.add_parser('illumina', help = 'Simulate haplotype sequences for three individuals')
    illumina.add_argument('fasta', type=str, help='Original FASTA sequence file.')
    illumina.add_argument('het', type=str, help = 'A FLOAT of heterozygosity level.')
    illumina.add_argument('dir', type = str, 
                          help = 'Output directory name.')

    pacbio = sub.add_parser('pacbio', help = 'Take simulated illumina data and generate pacbio reads')
    pacbio.add_argument('fastq', type = str, 
                        help = 'A FASTQ read file.')
    pacbio.add_argument('coverage', type = int, 
                        help = 'A FLOAT for the coverage in the pacbio read output. ')
    pacbio.add_argument('dir', type = str, 
                        help = 'Output directory name.')
    pacbio.add_argument('fasta', type = str, nargs = 6,
                        help = 'Original FASTA sequence files for each haplotype in trio. Please input them in the order of: MOM1 MOM2 DAD1 DAD2 CHILD1 CHILD2')

    pacbioccs = sub.add_parser('pacbioccs', help = 'Simulate haplotypes and corresponding pacbioccs reads')
    pacbioccs.add_argument('fasta', type=str,
                         help='Original FASTA haplotype sequence file.')
    pacbioccs.add_argument('het', type = str, help = "A FLOAT for heterozygosity level")
    pacbioccs.add_argument('dir', type = str, help = 'Output directory name.')

def main(sub, args):
    if sub == None:
        print(helpMessage)
    sim(sub, args)
