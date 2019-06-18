# WHdenovo
A cost-effective approach to diploid assembly for single samples and trios. It includes the following steps: construct sequence graph from illumina, align long reads to the graph and partition these long reads to two haplotypes.

Installation instructions: For starting, you need to compile some dependencies.

For SPAdes

```
cd trioasm/SPAdes-3.13.0
./spades_compile.sh
```

For whatshap
```
cd trioasm/whatshap
python setup.py build_ext -i
cd ../whatshap_trioasm
python setup.py build_ext -i
```
Other required dependencies are bfc, parallel and GraphAligner, which can be installed with ```conda indall```

For some python packags, we suggest manully copying them into your local python environment.
For example, if you are using Anaconda Python3.6:
```
cd trioasm/dependence
cp decorator.py ~/anaconda3/lib/python3.6/site-packages/
cp -r google ~/anaconda3/lib/python3.6/site-packages/
cp -r networkx ~/anaconda3/lib/python3.6/site-packages/
cp -r stream ~/anaconda3/lib/python3.6/site-packages/
cp -r vcf ~/anaconda3/lib/python3.6/site-packages/
cp -r xopen-0.3.5.dist-info ~/anaconda3/lib/python3.6/site-packages/
```


For simulating data
Illumina:
```
python src/simulate.py illumina N-bp.fasta <het> out/illumina/
```
And it will show you the file you'll need for pacbio data simulation and WHdenovo test
PacBio:
```
python src/simulate.py pacbio sample.fastq <coverage> out/pacbio/ <mom_hap1.fasta> <mom_hap2.fasta> <dad_hap1.fasta> <dad_hap2.fasta> <child_hap1.fasta> <child_hap2.fasta>
```

For running assembly
Trio case:
```
python src/assemble.py --illumina1 <illumina_dad_1.fq> --illumina2 <illumina_dad_2.fq> \
                       --pacbio <pacbio_mom.fasta> <pacbio_dad.fasta> <pacbio_child.fasta>
                       -p ped [-t <thread>] [-o out/path]
```

e.g.

```
python src/assemble.py --illumina1 16513_illumina_10k_1.5/child.het1.5.cov30_1.fq --illumina2 16513_illumina_10k_1.5/child.het1.5.cov30_2.fq \
                       --pacbio 16513_pacbio_10k_1.5_20/pacbio_mom.fasta 16513_pacbio_10k_1.5_20/pacbio_dad.fasta 16513_pacbio_10k_1.5_20/pacbio_child.fasta \
                       -p ped -t 24 -o test.simu
```

Individual case
```
python src/assemble.py --illumina1 <illumina_who_1.fq> --illumina2 <illumina_who_2.fq> \ 
                       --pacbio <pacbio_who.fasta> [-t <thread>] [-o out/path]
```
For assembling the genome from partitioned reads:
```
python set_dip_asm.py -f son.inputreads.fa -0 path/to/output/HP0.reads -1 path/to/output/HP1.reads --assemble -s 15k -t 40
```
Flye is the preset assembler in our pipeline. If you wish to use other assemblers with partitioned reads, just remove the ```--assemble``` argument.

For validating the partitioning of simulated the data.
```
python src/validate.py temp_pid_mmddhhmmss/bc1 <pacbio_who.fasta>
```
For validating the partitioning of real data when you have ground truth classification:
```
commands
```

We acknowledge the support of dependencies such as SPAdes, vg and GraphAligner.
