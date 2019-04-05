#!/bin/bash
#SOLiD test examples

art_SOLiD=../art_SOLiD
map2bed=../map2bed.pl
#art_SOLiD=../../bin/MacOS64/art_SOLiD 

# 1) singl-end 25bp reads simulation at 10X coverage 
$art_SOLiD -s testSeq.fa ./single_dat 25 10
echo "convert a map file to a UCSC BED file"
echo "$map2bed single_dat.bed single_dat.map"
$map2bed single_dat.bed single_dat.map

# 2) using user's read error profile for singl-end 75bp reads simulation at 20X coverage 
$art_SOLiD -s -p ../SOLiD_profiles/profile_pseudo ./testSeq.fa ./dat_userProfile 75 20
# 3) matepair 35bp (F3-R3) reads simulation at 20X coverage with DNA MEAN fragment size 2000bp and STD 50
$art_SOLiD -s testSeq.fa ./matepair_dat 35 20 2000 50
echo "convert two map files to a UCSC BED file"
echo "$map2bed maptepair.bed matepair_dat1.map matepair_dat2.map"
$map2bed matepair.bed matepair_dat_R3.map matepair_dat_F3.map

# 4) matepair reads simulation with a fixed random seed
$art_SOLiD -r 777 -s testSeq.fa ./matepair_fs1 50 10 1500 50
$art_SOLiD -r 777 -s testSeq.fa ./matepair_fs2 50 10 1500 50
echo "compare two simulation datasets"
diff matepair_fs1.sam matepair_fs2.sam

# 5) paired-end reads (75bp F3, 35bp F5) simulation with the MEAN fragment size 250 and STD 10 at 20X coverage
$art_SOLiD -s testSeq.fa ./paired_dat 75 35 50 250 10

# 6) amplicon sequencing with 25bp single-end reads at 100 reads per amplicon
$art_SOLiD -A s -s amplicon_reference.fa ./amp_single 25 100

# 7) amplicon sequencing with 50bp matepair reads at 80 read pairs per amplicon
$art_SOLiD -A m -s amplicon_reference.fa ./amp_matepair 50 80

# 8) amplicon sequencing with paired-end reads (35bp F3, 25bp F5 reads) at 50 pairs per amplicon
$art_SOLiD -A p -s amplicon_reference.fa ./amp_paired 35 25 50


