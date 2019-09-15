# WHdenovo
A cost-effective approach to diploid assembly for single samples and trios. It includes the following steps: construct sequence graph from illumina, align long reads to the graph and partition these long reads to two haplotypes.

### Installation

Conda installation
```
conda install -c bioconda whdenovo
```

Developer installation
```
wget <release>.tar.gz
tar -xvzf <release>.tar.gz
cd WHdenovo
python setup.py install
```
If 

### Data Simulation

For simulating Illumina data:
```
whdenovo simulate illumina N-bp.fasta <het> out/illumina/
```
And it will show you the file you'll need for pacbio data simulation and WHdenovo test

For simulating PacBio data:
```
whdenovo simulate pacbio sample.fastq <coverage> out/pacbio/ \
                       <mom_hap1.fasta> <mom_hap2.fasta> \
                       <dad_hap1.fasta> <dad_hap2.fasta> \
                       <child_hap1.fasta> <child_hap2.fasta>
```
### Assembly

Trio case:
```
whdenovo partition --illumina1 <illumina_child_1.fq> --illumina2 <illumina_child_2.fq> \
                       --pacbio <pacbio_mom.fasta> <pacbio_dad.fasta> <pacbio_child.fasta>
                       -p ped [-t <thread>] [-o out/path] [--lowc INT] [--highc INT]
```

e.g.

```
whdenovo partition --illumina1 16513_illumina_10k_1.5/child.het1.5.cov30_1.fq --illumina2 16513_illumina_10k_1.5/child.het1.5.cov30_2.fq \
                       --pacbio 16513_pacbio_10k_1.5_20/pacbio_mom.fasta 16513_pacbio_10k_1.5_20/pacbio_dad.fasta 16513_pacbio_10k_1.5_20/pacbio_child.fasta \
                       -p ped -t 24 -o test.simu --lowc 5 --highc 60
```

Individual case
```
whdenovo partition --illumina1 <illumina_who_1.fq> --illumina2 <illumina_who_2.fq> \ 
                       --pacbio <pacbio_who.fasta> [-t <thread>] [-o out/path]
```
For assembling the genome from partitioned reads given by ```partition.py```:

```
conda activate flye
whdenovo assemble -f son.inputreads.fa \
                      -0 path/to/output/HP0.reads -1 path/to/output/HP1.reads \
                      --assemble -s 15k -t 40 <--pacbio|--nano>
```
If you wish to use other assemblers with partitioned reads, just remove the ```--assemble``` argument and you may not need to activate the virtual environment.

### Result Validation

For validating the partitioning of simulated the data.
```
whdenovo validate -p out/path -f <pacbio_who.fasta>
```
For validating the partitioning of real data when you have ground truth classification:
```
whdenovo validate -p out/path -t <tagged.reads.txt>
```
An example for tagged.reads.txt is at test/haplotagged.reads.txt, which should include the HP tag and PS tag from haplotagged BAM file.
***
We acknowledge the support of dependencies such as [bfc](https://github.com/lh3/bfc), [SPAdes](http://cab.spbu.ru/software/spades/), [vg](https://github.com/vgteam/vg) and [GraphAligner](https://github.com/maickrau/GraphAligner).
