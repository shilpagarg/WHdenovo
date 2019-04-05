#!/bin/bash
#illumina test examples
art=../art_illumina
#art=../../bin/MacOS64/art_illumina
#art=../../bin/Linux64/art_illumina

# 1) simulation of single-end reads of 35bp with 10X using the built-in combined quality profile, and without Indels
echo 11111
$art -ss GA1 -i ./testSeq.fa -o ./single_end_com -l 35 -f 10 -sam -k 0
#convert an aln file to a bed file
../aln2bed.pl single_end_com.bed single_end_com.aln
echo 22222
# 2) simulation of single-end reads of 50bp with 10X using the built-in seperated quality profiles for A, C, G, and T 
$art -ss MinS -i ./testSeq.fa -o ./single_end_sep -l 50 -f 10 -sp -sam
#convert an aln file to a bed file
../aln2bed.pl single_end_sep.bed single_end_sep.aln
echo 33333
# 3) simulation of paired-end reads of 150bp with the mean fragment size 500 and standard deviation 10
#    using the built-in combined read quality profiles
$art -ss HS25 -i ./testSeq.fa -o ./paired_end_com -l 150 -f 10 -p -m 500 -s 10 -sam
#convert both aln files to a bed file
../aln2bed.pl paired_end_com.bed paired_end_com1.aln paired_end_com2.aln
echo 44444
# 4) simulation of paired-end reads of 100bp with the mean fragment size 500 and standard deviation 10
#   using the built-in seperated quality profiles for A, C, G, and T 
$art -ss HS20 -i ./testSeq.fa -o ./paired_end_sep -l 100 -f 10 -p -m 500 -s 10 -sp -sam
#convert both aln files to a bed file
../aln2bed.pl paired_end_sep.bed paired_end_sep1.aln paired_end_sep2.aln
echo 55555
# 5) simulation of mate-pair reads of 100bp with the mean fragment size 2500 and standard deviation 50
#    using the built-in combined read quality profiles
$art -ss HS25 -i ./testSeq.fa -o ./matepair_com -l 100 -f 10 -p -m 2500 -s 50 -sam
#convert both aln files to a bed file
../aln2bed.pl matepair_com.bed matepair_com1.aln matepair_com2.aln
echo 66666
# 6) amplicaton read simulation: generate two 100bp single-end reads from 5' end for each amplicon reference
$art -ss HSXt -i ./amplicon_reference.fa  -amp  -o ./amp_5_end_com -l 100 -f 2 -sam
#convert both aln files to a bed file
../aln2bed.pl amp_5_end_com.bed  amp_5_end_com.aln
echo 77777
# 7) amplicaton read simulation: generate one 100bp paired-end reads from both ends for each amplicon reference
$art -ss MSv1 -i ./amplicon_reference.fa  -amp  -o ./amp_pair -p -l 100 -f 1 -sam
#convert both aln files to a bed file
../aln2bed.pl amp_pair.bed amp_pair1.aln amp_pair2.aln
echo 88888
# 8) amplicaton read simulation: generate one 100bp matepair reads from both ends for each amplicon reference
$art -ss MSv3 -i ./amplicon_reference.fa  -amp  -o ./amp_matepair -mp -l 100 -f 1 -sam
#convert both aln files to a bed file
../aln2bed.pl amp_matepair.bed amp_matepair1.aln amp_matepair2.aln
echo 99999
# 9) generate two identical simulation datasets with a fixed random seed
$art -ss HSXn -i ./testSeq.fa -rs 777 -o ./paired_end_com_f1 -l 100 -f 10 -p -m 500 -s 10 -sam
$art -ss HSXn -i ./testSeq.fa -rs 777 -o ./paired_end_com_f2 -l 100 -f 10 -p -m 500 -s 10 -sam
echo 00000
# 10) reduce the sequencing error rate to one 10th of the default profile for a paired-end read simulation  
$art -ss HS20 -i ./testSeq.fa -qs 10 -qs2 10 -o ./paired_end_com_f1 -l 100 -f 10 -p -m 500 -s 10 -sam
echo 1111111111
# 11) turn off the masking of 'N' genomic regions  
$art -ss NS50 -nf 0 -i ./testSeq.fa -o ./paired_nomask -l 75 -f 10 -p -m 500 -s 10 -sam
