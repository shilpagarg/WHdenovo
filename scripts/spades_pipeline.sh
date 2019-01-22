
# combine end reads for both ends into same files 
cat trioyeast/DBVPG6044.R1.fastq trioyeast/DBVPG6765.R1.fastq trioyeast/Y12.R1.fastq trioyeast/YPS128.R1.fastq> illumina_1_total.fq
cat trioyeast/DBVPG6044.R2.fastq trioyeast/DBVPG6765.R2.fastq trioyeast/Y12.R2.fastq trioyeast/YPS128.R2.fastq> illumina_2_total.fq

python2 /Users/IsaacSebenius/try/spades.py -t 16 -1 illumina_1_total.fq -2 illumina_2_total.fq -k 21 --only-assembler -o spades

grep -v '^P' spades/assembly_graph.gfa | awk -F'\t' '{ if ($2 != $4) print $0}' | vg view --gfa-in - --vg | vg view -g - | awk -F'\t' '{ if ($2 !=$4) print $0}' > graph.gfa

# #get rid of nodes with degree 0 or over 50
python2 /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/printnodedegrees_gfa.py graph.gfa | awk -F' ' '{ if($2 > 50 || $2==0) printf "%s\n", $1 }' > wrongnodes.txt
python /Users/IsaacSebenius/Dropbox/HARVARD/GargLab/Yeast/remove_wrongnodes.py wrongnodes.txt graph.gfa