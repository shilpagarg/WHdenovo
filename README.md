# WHdenovo
A cost-effective approach to diploid assembly for single samples and trios. It includes the following steps: construct sequence graph from illumina, align long reads to the graph and partition these long reads to two haplotypes.

### Installation

For cloning this repository, you will need the help of [Git Large File Storage](https://github.com/git-lfs/git-lfs), because of some usable pre-compiled dependencies.

This can be installed in various ways:
```
sudo yum install -y git-lfs
```
```
sudo apt-get install git-lfs
```
Then clone this repository with:
```
git lfs clone https://github.com/shilpagarg/WHdenovo.git
```

For starting, you need to compile some dependencies.

For SPAdes

```
cd trioasm/SPAdes-3.13.0
./spades_compile.sh
```

For WhatsHap
```
cd trioasm/whatshap
python setup.py build_ext -i
cd ../whatshap_trioasm
python setup.py build_ext -i
```
Other required dependencies are [bfc](https://github.com/lh3/bfc), [parallel](https://www.gnu.org/software/parallel/) and [GraphAligner](https://github.com/maickrau/GraphAligner), which can be installed with ```conda install```. We adopt [Flye](https://github.com/fenderglass/Flye) as preset for the final assembly. Flye requires Python 2.x environment while our main pipeline requires python 3.x. Thus, before installing Flye, we strongly suggest preparing a virtual environment for it ```conda create --name flye```, and then ```conda install flye```.

For some python packages, we suggest manully copying them into your local python environment.
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

### Data Simulation

For simulating Illumina data:
```
python src/simulate.py illumina N-bp.fasta <het> out/illumina/
```
And it will show you the file you'll need for pacbio data simulation and WHdenovo test

For simulating PacBio data:
```
python src/simulate.py pacbio sample.fastq <coverage> out/pacbio/ \
                       <mom_hap1.fasta> <mom_hap2.fasta> \
                       <dad_hap1.fasta> <dad_hap2.fasta> \
                       <child_hap1.fasta> <child_hap2.fasta>
```
### Assembly

Trio case:
```
python src/partition.py --illumina1 <illumina_child_1.fq> --illumina2 <illumina_child_2.fq> \
                       --pacbio <pacbio_mom.fasta> <pacbio_dad.fasta> <pacbio_child.fasta>
                       -p ped [-t <thread>] [-o out/path]
```

e.g.

```
python src/partition.py --illumina1 16513_illumina_10k_1.5/child.het1.5.cov30_1.fq --illumina2 16513_illumina_10k_1.5/child.het1.5.cov30_2.fq \
                       --pacbio 16513_pacbio_10k_1.5_20/pacbio_mom.fasta 16513_pacbio_10k_1.5_20/pacbio_dad.fasta 16513_pacbio_10k_1.5_20/pacbio_child.fasta \
                       -p ped -t 24 -o test.simu
```

Individual case
```
python src/partition.py --illumina1 <illumina_who_1.fq> --illumina2 <illumina_who_2.fq> \ 
                       --pacbio <pacbio_who.fasta> [-t <thread>] [-o out/path]
```
For assembling the genome from partitioned reads given by ```partition.py```:
```
conda activate flye
python set_dip_asm.py -f son.inputreads.fa \
                      -0 path/to/output/HP0.reads -1 path/to/output/HP1.reads \
                      --assemble -s 15k -t 40
```
If you wish to use other assemblers with partitioned reads, just remove the ```--assemble``` argument and you may not need to activate the virtual environment.

### Result Validation

For validating the partitioning of simulated the data.
```
python src/validate.py out/path/bc1 <pacbio_who.fasta>
```
For validating the partitioning of real data when you have ground truth classification:
```
python src/validate_real.py out/path/bc1 <tagged.reads.txt>
```
An example for tagged.reads.txt is at test/haplotagged.reads.txt, which should include the HP tag and PS tag from haplotagged BAM file.
***
We acknowledge the support of dependencies such as [bfc](https://github.com/lh3/bfc), [SPAdes](http://cab.spbu.ru/software/spades/), [vg](https://github.com/vgteam/vg) and [GraphAligner](https://github.com/maickrau/GraphAligner).
