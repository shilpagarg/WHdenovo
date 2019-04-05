#USAGE: selected_partitioned_reads all_partitioned_reads all_pacbio_reads.fq true_hap1.fa true_hap2.fa genomeSize 
echo 
echo PARTITIONING STATS FOR SELECTED READS
echo ======================================
#print partitioning accuracy for selected

python /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/compute_read_partitioning_test.py $1
python /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/compute_read_partitioning_test.py $1 > selected_partitioning_accuracy.txt

echo 
echo PARTITIONING STATS FOR ALL PARTITIONED READS
echo ============================================

#get partitioning accuracy for all partitioned reads
python /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/compute_read_partitioning_all_test.py $2 
python /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/compute_read_partitioning_all_test.py $2 > all_partitioning_accuracy.txt

echo 
echo WORKING ON GETTING ASSMEMLY...
echo ==============================

awk '{if ($3 == 0) print $1}' $2 | grep 'momHAP' > mom_part1.txt
awk '{if ($3 == 0) print $1}' $2 | grep 'dadHAP' > dad_part1.txt
awk '{if ($3 == 0) print $1}' $2 | grep 'childHAP' > child_part1.txt
awk '{if ($3 == 1) print $1}' $2 | grep 'momHAP' > mom_part2.txt
awk '{if ($3 == 1) print $1}' $2 | grep 'dadHAP' > dad_part2.txt
awk '{if ($3 == 1) print $1}' $2 | grep 'childHAP' > child_part2.txt

#print and store partitioning accuracy 
python /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/Yeast/get_unpartitioned.py $2 $3
cat child_part1.txt mom_part1.txt unpartitioned.txt > part_1_read_names.txt
cat child_part2.txt dad_part1.txt unpartitioned.txt > part_2_read_names.txt

/Users/IsaacSebenius/seqtk/seqtk subseq $3 part_1_read_names.txt > part_1_reads.fq
/Users/IsaacSebenius/seqtk/seqtk subseq $3 part_2_read_names.txt > part_2_reads.fq

#run canu to do assembly
/Users/IsaacSebenius/canu/canu/Darwin-amd64/bin/canu -p asm -d hap1 genomeSize=$6 stopOnLowCoverage=10 -pacbio-raw part_1_reads.fq

echo 
echo FIRST HAP DONE
echo ==============
/Users/IsaacSebenius/canu/canu/Darwin-amd64/bin/canu -p asm -d hap2 genomeSize=$6 stopOnLowCoverage=10 -pacbio-raw part_2_reads.fq

echo
echo SECOND HAP DONE
echo ===============

echo 
echo GENERATING ASSEMBLY STATS
echo =========================

#run countfasta
perl /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/countFasta.pl /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/Yeast/new_chrom1_test/hap1/asm.contigs.fasta > hap1_count_fasta_stats.txt

#¡¡¡¡ NOTE: need to adjust the following commands based on what actually is the true haplotype !!!! - for this example I'm not sure if this is right, this just shows the pipeline goes ok. 

#run dnadiff
/Users/IsaacSebenius/anaconda3/opt/mummer-3.23/dnadiff $4 /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/Yeast/new_chrom1_test/hap1/asm.contigs.fasta
mv out.report hap1_dna_diff_output.txt

perl /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/countFasta.pl /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/Yeast/new_chrom1_test/hap2/asm.contigs.fasta > hap2_count_fasta_stats.txt
/Users/IsaacSebenius/anaconda3/opt/mummer-3.23/dnadiff $5 /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/Yeast/new_chrom1_test/hap2/asm.contigs.fasta
mv out.report hap2_dna_diff_output.txt

echo 
echo GOODBYE!
echo ========

#child2 is dad1, child1 is mom1