# WHdenovo
Another method for diploid assembly, which works only with PacBioCCS data.

```
usage: haplotype.py [-h] [-p PAF] [-o GFA] [-t INT] [--het FLOAT]
                          [-c INT]
                          FASTA
positional arguments:
  FASTA                 Input all read records in one FASTA file.

optional arguments:
  -h, --help            show this help message and exit
  -p PAF, --paf PAF     PAF file, should be consistent with the input FASTA
                        file; if not given, we will run minimap2 on the input
                        fasta and generate one.
  -o GFA, --output GFA  Output to FILE, default to stdout.
  -t INT, --threads INT
                        Maximum number of threads to use. [4]
  --het FLOAT           Heterozygosity rate, set as a filter for abnormal
                        alignment. [0.1]
  -c INT, --avg-coverage INT
                        average coverage, set as a filter for abnormal
                        alignment. [40]
```
In the output GFA graph, we include an extra column in the "L" lines (for edges) with one of the color tags "GREEN"/"RED"/"DASH", which stand for belonging to the same haplotype, different haplotype, and don't have enough phasing information but overlap, repectively.

An example run on our test data:
```
python haplotype.py -t 40 test.fasta > test.gfa
```
For visualization of the graph, you may extract only the GREEN edges to have a clear view of the partitioning.
```
grep -E '^S|GREEN' test.gfa > test.green.gfa
```
