#!/bin/bash
#454 test examples
art_454=../art_454
#art_454=../../bin/MacOS64/art_454
#art_454=../../bin/Linux64/art_454

# 1) simulation of 454 single-end reads at 10X coverage using the built-in 454 FLX read  profile
$art_454 -a -s ./testSeq.fa ./single_454_flx 10
#   convert an aln file to a bed file
../aln2bed.pl single_454_flx.bed single_454_flx.aln

# 2) simulation of 454 paried-end reads at 5X coverage with the mean fragment size 500 
#    and standard deviation 20 using the built-in 454 FLX read  profile
$art_454 -a -s ./testSeq.fa ./paired_454_flx 5 500 20
#    convert both aln files to a bed file
../aln2bed.pl paired_454_flx.bed paired_454_flx1.aln paired_454_flx2.aln

# 3) simulation of 454 paried-end reads at 6X coverage with the mean fragment size 500 
#    and standard deviation 20 using the built-in 454 FLX Titanium read  profile
$art_454 -a -s -t ./testSeq.fa ./paired_454_flxTitan 6 500 20
#   convert both aln files to a bed file
../aln2bed.pl paired_454_flxTitan.bed paired_454_flxTitan1.aln paired_454_flxTitan2.aln 

# 4) 5' end amplicon sequencing with 10 single-end reads per amplicon   
$art_454 -a -s -A ./amplicon_reference.fa ./amp_single_454 10

# 5) both 5' and 3' ends amplicon sequencing with 5 read pairs per amplicon  
$art_454 -a -s -B ./amplicon_reference.fa ./amp_paired_454 5 

# 6) used a fixed random seed to generate identical datasets   
$art_454 -a -s -r 777 ./testSeq.fa ./single_454_t1 10
$art_454 -a -s -r 777 ./testSeq.fa ./single_454_t2 10
# compare difference
echo "compare difference of two simulation datasets"
diff  single_454_t1.sam  single_454_t2.sam

