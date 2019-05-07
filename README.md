# WHdenovo
A cost-effective approach to diploid assembly for single samples and trios. It includes the following steps: construct sequence graph from illumina, align long reads to the graph and partition these long reads to two haplotypes.

Installation instructions: For starting, you need to compile some dependencies.

for SPAdes

```
cd trioasm/SPAdes-3.13.0
./spades_compile.sh
```

for whatshap
```
cd trioasm/whatshap
python setup.py build_ext -i
cd ../whatshap_trioasm
python setup.py build_ext -i
```

And for some python packags, you need to manully copy them into your local python environment.
For example, if you are using Anaconda Python3.6:
```
cd trioasm/dependence
cp decorator.py ~/anaconda3/lib/python3.6/site-packages/
cp -r google ~/anaconda3/lib/python3.6/site-packages/
cp -r networkx ~/anaconda3/lib/python3.6/site-packages/
cp -r stream ~/anaconda3/lib/python3.6/site-packages/
cp -r vcf ~/anaconda3/lib/python3.6/site-packages/
cp -r xopen-0.3.5.dist-info
``` 


For simulating data
Illumina:
```
python src/simulate.py illumina N-bp.fasta <het> out/illumina/
```
And it will show you the file you'll need
```
python src/simulate.py pacbio sample.fastq <coverage> out/pacbio/ <mom_hap1.fasta> <mom_hap2.fasta> <dad_hap1.fasta> <dad_hap2.fasta> <child_hap1.fasta> <child_hap2.fasta>
```
You can directly take the illumina data the program told you to use when simulating illumina data.


For running assembly
Trio case. Illumina file names and paths are provided when you simulate.

```
python src/assemble.py --illumina1 <illumina_mom_1.fq> <illumina_dad_1.fq> <illumina_dad_1.fq> \
                       --illumina2  <illumina_mom_2.fq> <illumina_dad_2.fq> <illumina_dad_2.fq> \
                       --pacbio <pacbio_mom.fasta> <pacbio_dad.fasta> <pacbio_child.fasta>  -p ped \ 
                       -t <thread> 
```

e.g.

```
python src/assemble.py --illumina1 16513_illumina_10k_1.5/mom.het1.5.cov30_1.fq 16513_illumina_10k_1.5/dad.het1.5.cov30_1.fq 16513_illumina_10k_1.5/child.het1.5.cov30_1.fq \
                       --illumina2 16513_illumina_10k_1.5/mom.het1.5.cov30_2.fq 16513_illumina_10k_1.5/dad.het1.5.cov30_2.fq 16513_illumina_10k_1.5/child.het1.5.cov30_2.fq \
                       --pacbio 16513_pacbio_10k_1.5_20/pacbio_mom.fasta 16513_pacbio_10k_1.5_20/pacbio_dad.fasta 16513_pacbio_10k_1.5_20/pacbio_child.fasta \
                       -p ped -t 24
```


Individual case
```
python src/assemble.py --illumina1 <illumina_who_1.fq> --illumina2 <illumina_who_2.fq> \ 
                       --pacbio <pacbio_who.fasta> -t <thread>
```

For validation
First you need to know which temp directory you just output
```
python src/validate.py temp_pid_mmddhhmmss/bc1 <pacbio_who.fasta>

```
