ART_ILLUMINA  README (updated on 06/6/2016) Weichun Huang <whduke@gmail.com>
ART_Illumina (2008-2016), Q Version 2.5.8 (Jun 6, 2016)                      

DESCRIPTION

	ART (art_Illumina Q version) is a simulation program to generate sequence read data of Illumina
       	sequencers. ART generates reads according to the empirical read quality profile summarized from
       	large real read data. ART has been using for testing or benchmarking a variety of methods or tools
       	for next-generation sequencing data analysis, including read alignment, de novo assembly, detection
        of SNP, CNV, or other structure variation.

	art_Illumina can generate single-end, paired-end,  mate-pair reads, and amplicon sequencing simulation
       	of Illumina sequencing platform. Its outputs include FASTQ read, ALN and/or SAM alignment files. ALN
       	files can also be converted to UCSC BED files by using aln2bed.pl program included. 

	art_Illumina comes with the tool art_profiler_illumina that can generate quality profiles from Illumina
        sequencing data in the fastq format. The tool is in the folder ART_profiler_illumina. Please see README
       	in the folder for the details and usage.

COMPILATION AND INSTALLATION

	PREREQUISITES: 

		1) GNU g++ 4.0 or above (http://gcc.gnu.org/install) 
		2) GNU gsl library (http://www.gnu.org/s/gsl/)

	COMPILATION & INSTALLATION
		1) add gsl library installation directory to compiler's search path

	       	If your GSL library is not installed in system standard search path, your GSL installation 
		directory need to be added to gcc/g++ complier FLAGS. You can find your gsl installation path by
	       	running the command "gsl-config —prefix” in your terminal window. gsl is typically installed under
		path /usr/local or /opt/local in MacOS X. You can add your gsl installation path to g++ searching
	       	path by running the following command:

		export CFLAGS="$CFLAGS -I/opt/local/include" CPPFLAGS="$CPPFLAGS -I/opt/local/include" LDFLAGS="$LDFLAGS -L/opt/local/lib"
		export CFLAGS="$CFLAGS -I/usr/local/include" CPPFLAGS="$CPPFLAGS -I/usr/local/include" LDFLAGS="$LDFLAGS -L/usr/local/lib"

		2) compile and install 

		./configure --prefix=$HOME
	       	make
	       	make install

EXAMPLES

	In the "examples" subdirectory, the shell script "run_test_examples_illumina.sh" gives four examples of using
       	ART for read simulation.  To test these four examples, just run the script "run_test_examples_illumina.sh"

USAGE
	RECOMMENDED USAGES (specifying a sequencing system)
		art_illumina [options] -ss <sequencing_system> -sam -i <seq_ref_file> -l <read_length> -f <fold_coverage> -o <outfile_prefix>
		art_illumina [options] -ss <sequencing_system> -sam -i <seq_ref_file> -l <read_length> -c <num_reads_per_sequence> -o <outfile_prefix>
		art_illumina [options] -ss <sequencing_system> -sam -i <seq_ref_file> -l <read_length> -f <fold_coverage> -m <mean_fragsize> -s <std_fragsize> -o <outfile_prefix>
		art_illumina [options] -ss <sequencing_system> -sam -i <seq_ref_file> -l <read_length> -c <num_reads_per_sequence> -m <mean_fragsize> -s <std_fragsize> -o <outfile_prefix>

	OTHER USAGES
		art_illumina [options] -sam -i <seq_ref_file> -l <read_length> -f <fold_coverage> -o <outfile_prefix>
		art_illumina [options] -sam -i <seq_ref_file> -l <read_length> -c <total_num_reads> -o <outfile_prefix>
		art_illumina [options] -sam -i <seq_ref_file> -l <read_length> -f <fold_coverage> -m <mean_fragsize> -s <std_fragsize> -o <outfile_prefix>
		art_illumina [options] -sam -i <seq_ref_file> -l <read_length> -c <total_num_reads> -m <mean_fragsize> -s <std_fragsize> -o <outfile_prefix>
	
	===== PARAMETERS =====
	    -1   --qprof1   the first-read quality profile
	    -2   --qprof2   the second-read quality profile
	    -amp --amplicon amplicon sequencing simulation
	    -c   --rcount   number of reads/read pairs to be generated per sequence(not be used together with -f/--fcov)
	    -d   --id       the prefix identification tag for read ID
	    -ef  --errfree  indicate to generate the zero sequencing errors SAM file as well the regular one
	                    NOTE: the reads in the zero-error SAM file have the same alignment positions
	                    as those in the regular SAM file, but have no sequencing errors
	    -f   --fcov     the fold of read coverage to be simulated or number of reads/read pairs generated for each amplicon
	    -h   --help     print out usage information
	    -i   --in       the filename of input DNA/RNA reference
	    -ir  --insRate  the first-read insertion rate (default: 0.00009)
	    -ir2 --insRate2 the second-read insertion rate (default: 0.00015)
	    -dr  --delRate  the first-read deletion rate (default:  0.00011)
	    -dr2 --delRate2 the second-read deletion rate (default: 0.00023)
	    -k   --maxIndel the maximum total number of insertion and deletion per read (default: up to read length)
	    -l   --len      the length of reads to be simulated
	    -m   --mflen    the mean size of DNA/RNA fragments for paired-end simulations
	    -mp  --matepair indicate a mate-pair read simulation
	    -M  --cigarM    indicate to use CIGAR 'M' instead of '=/X' for alignment match/mismatch
	    -nf  --maskN    the cutoff frequency of 'N' in a window size of the read length for masking genomic regions
	                    NOTE: default: '-nf 1' to mask all regions with 'N'. Use '-nf 0' to turn off masking
	    -na  --noALN    do not output ALN alignment file
	    -o   --out      the prefix of output filename
	    -p   --paired   indicate a paired-end read simulation or to generate reads from both ends of amplicons
	                    NOTE: art will automatically switch to a mate-pair simulation if the given mean fragment size >= 2000
	    -q   --quiet    turn off end of run summary
	    -qL  --minQ     the minimum base quality score
	    -qU  --maxQ     the maxiumum base quality score
	    -qs  --qShift   the amount to shift every first-read quality score by 
	    -qs2 --qShift2  the amount to shift every second-read quality score by
	                    NOTE: For -qs/-qs2 option, a positive number will shift up quality scores (the max is 93) 
	                    that reduce substitution sequencing errors and a negative number will shift down 
	                    quality scores that increase sequencing errors. If shifting scores by x, the error
	                    rate will be 1/(10^(x/10)) of the default profile.
	    -rs  --rndSeed  the seed for random number generator (default: system time in second)
	                    NOTE: using a fixed seed to generate two identical datasets from different runs
	    -s   --sdev     the standard deviation of DNA/RNA fragment size for paired-end simulations.
	    -sam --samout   indicate to generate SAM alignment file
	    -sp  --sepProf  indicate to use separate quality profiles for different bases (ATGC)
	    -ss  --seqSys   The name of Illumina sequencing system of the built-in profile used for simulation

	    NOTE: all built-in sequencing system ID names are:
		GA1 - GenomeAnalyzer I (36bp,44bp)
		GA2 - GenomeAnalyzer II (50bp, 75bp) 
		HS10 - HiSeq 1000 (100bp)
		HS20 - HiSeq 2000 (100bp)
		HS25 - HiSeq 2500 (125bp, 150bp)
		HSXn - HiSeqX PCR free (150bp)
		HSXt - HiSeqX TruSeq (150bp)
		MinS - MiniSeq TruSeq (50bp)
		MSv1 - MiSeq v1 (250bp)
		MSv3 - MiSeq v3 (250bp)
		NS50 - NextSeq500 v2 (75bp)

	  ===== NOTES =====
	  
	  * ART by default selects a built-in quality score profile according to the read length specified for the run.
	  
	  * For single-end simulation, ART requires input sequence file, outputfile prefix, read length, and read count/fold coverage.
	  
	  * For paired-end simulation (except for amplicon sequencing), ART also requires the parameter values of
	    the mean and standard deviation of DNA/RNA fragment lengths
	  
	  ===== EXAMPLES =====
	  
	   1) single-end read simulation
	   	art_illumina -ss HS25 -sam -i reference.fa -l 150 -f 10 -o single_dat
	  
	   2) paired-end read simulation
	         art_illumina -ss HS25 -sam -i reference.fa -p -l 150 -f 20 -m 200 -s 10 -o paired_dat
	  
	   3) mate-pair read simulation
	         art_illumina -ss HS10 -sam -i reference.fa -mp -l 100 -f 20 -m 2500 -s 50 -o matepair_dat
	  
	   4) amplicon sequencing simulation with 5' end single-end reads 
	   	art_illumina -ss GA2 -amp -sam -na -i amp_reference.fa -l 50 -f 10 -o amplicon_5end_dat
	  
	   5) amplicon sequencing simulation with paired-end reads
	         art_illumina -ss GA2 -amp -p -sam -na -i amp_reference.fa -l 50 -f 10 -o amplicon_pair_dat
	  
	   6) amplicon sequencing simulation with matepair reads
	         art_illumina -ss MSv1 -amp -mp -sam -na -i amp_reference.fa -l 150 -f 10 -o amplicon_mate_dat
	  
	   7) generate an extra SAM file with zero-sequencing errors for a paired-end read simulation
	         art_illumina -ss HSXn -ef -i reference.fa -p -l 150 -f 20 -m 200 -s 10 -o paired_twosam_dat
	  
	   8) reduce the substitution error rate to one 10th of the default profile
	         art_illumina -i reference.fa -qs 10 -qs2 10 -l 50 -f 10 -p -m 500 -s 10 -sam -o reduce_error
	  
	   9) turn off the masking of genomic regions with unknown nucleotides 'N'
	         art_illumina -ss HS20 -nf 0  -sam -i reference.fa -p -l 100 -f 20 -m 200 -s 10 -o paired_nomask
	  
	   10) masking genomic regions with >=5 'N's within the read length 50
	         art_illumina -ss HSXt -nf 5 -sam -i reference.fa -p -l 150 -f 20 -m 200 -s 10 -o paired_maskN5
	  
READ QUALITY PROFILE

	TOOL FOR CREATING A NEW QUALITY PROFILE

        art_profiler_illumina in the folder ART_profiler_illumina can generate quality profiles from Illumina
	fastq data. Please see README in the folder for the usage.

	FORMAT 
	A valid quality score profile is tab-delimited and has no specific header line. Headers can be included
       	if each line of extraneous information begins with a number sign(#). Each line of actual quality profile
       	information must begin with an identifier that indicates where the data comes from for the remainder of
       	the line. Identifiers include:

		.	The identifier for the combination of all base information.
		
		A	The identifier for quality scores associated with A calls only.

		T	The identifier for quality scores associated with T calls only.

		G	The identifier for quality scores associated with G calls only.

		C	The identifier for quality scores associated with C calls only.


	Following the identifier on a given line must be the position number indicating where the rest of the
       	information on that line applies within a given fragment. The data must be arrayed in pairs such that each
       	line in the pair has the same identifier and position number.  The first line in a pair is a list of quality
       	scores in ascending order and the second line are the corresponding cumulative frequencies of the quality scores. 


	EXAMPLE:


	.       0       3       6       7       8       9       10      11      12	13      14      15
	.       0       39375   355755  395136  415685  1227131 1338634 1522001	1851208 2165909 2436839 2608538
	...
	.       35      4       5       6       7       8       9       10      11	12      13      14
	.       35      434262  1341151 1725690 2979293 3478620 3592624 3672807	3873754 4096922 4983957 6111261
	A       0       3       6       7       8       9       10      11      12	13      14      15
	A       0       560     99637   111485  119727  389899  412150  458066  572958	665153  745916  793532

	BUILT-IN PROFILES 

	Under subfolder "Illumina_Profiles" are the read profiles of Illumina GAII sequencers including all
	ART's built-in profiles as indicated below. 
      
	1) Raw quality profiles
		GA1 36bp reads
			Emp36R1.txt
	       		Emp36R2.txt

		GA1 44bp reads
			Emp44R1.txt
		       	Emp44R2.txt

		GA2 50bp reads
	       		Emp50R1.txt
		       	Emp50R2.txt

		GA2 75bp reads
	       		Emp75R1.txt
		       	Emp75R2.txt

		MiSeq 250bp reads (250bp reads)
       			EmpMiSeq250R1.txt 
       			EmpMiSeq250R2.txt 

		HiSeq1000 100bp reads 
	       		Emp100R1.txt
		       	Emp100R2.txt

		HiSeq2000 100bp reads (the new default profile for 100bp reads)
		       	HiSeq2kL100R1.txt
		     	HiSeq2kL100R2.txt

		HiSeq2500 125bp reads 
			HiSeq2500L125R1.txt
		       	HiSeq2500L125R2.txt

		HiSeq2500 150bp reads 
		       	HiSeq2500L150R1.txt
		       	HiSeq2500L150R2.txt

		HiSeqX PCR free (150bp)
			HiSeqXPCRfreeL150R1.txt
			HiSeqXPCRfreeL150R2.txt

		HiSeqX TruSeq (150bp)
			HiSeqXtruSeqL150R1.txt
			HiSeqXtruSeqL150R2.txt

		MiniSeq TruSeq (50bp)
			MiniSeqTruSeqL50.txt

		MiSeq v3 (250bp)
			MiSeqv3L250R1.txt
			MiSeqv3L250R2.txt

		NextSeq500 v2 (75bp)
			NextSeq500v2L75R1.txt
			NextSeq500v2L75R2.txt

	2) Recalibrated quality profiles (all these are ART's built-in profiles) 

		36bp reads
	       		1st: EmpR36R1
		       	2nd: EmpR36R2

		44bp reads
	       		1st: EmpR44R1
		       	2nd: EmpR44R2

	    	50bp reads
	       		1st: EmpR50R1
			2nd: EmpR50R2

	    	75 reads
			1st: EmpR75R1
			2nd: EmpR75R2

OUTPUT FILES

	*.fq   - FASTQ read data files. For paired‐read simulation, *1.fq contains data of the first reads, and *2.fq for the second reads.
	*.aln  - read alignment files. For paired‐read simulation, *1.aln has read alignments for the first reads and *2.aln for the second reads.

	FASTQ FILE format 
		A FASTQ file contains both sequence bases and quality scores of sequencing reads and is in the following format:  
			@read_id 
			sequence_read 
			+ 
			base_quality_scores 
	
		A base quality score is coded by the ASCII code of a single character, where the quality score is equal to ASCII code of the
       		character minus 33.    

		Example: 
		@refid-4028550-1 
		caacgccactcagcaatgatcggtttattcacgat...
		+ 
		????????????7?????<??>??=&?<<?-<?0?...

	ALN format 
		An ALN file has a Header and main Body parts. The header part includes the command used to generate this file and reference
	       	sequence id and length. The header @CM tag for command line, and @SQ for reference sequence.  A header always starts with 
		"##ART" and ends with  "##Header End".

		HEADER EXAMPLE

		##ART_Illumina  read_length     35
		@CM     ../art_illumina -i ./testSeq.fa -o ./single_end_com -l 35 -f 10 -sam -rs 177
		@SQ     seq1    7207
	       	@SQ     seq2    3056
		##Header End

		The body part contains all alignments in the following format 

		>ref_seq_id	read_id	aln_start_pos	ref_seq_strand
	       	ref_seq_aligned
	       	read_seq_aligned 
	
		aln_start_pos is the alignment start position of reference sequence. aln_start_pos is always relative to the strand of reference
       		sequence. That is, aln_start_pos 10 in the plus (+) strand is different from aln_start_pos 10 in the minus (‐) stand.  
	
		ref_seq_aligned is the aligned region of reference sequence, which can be from plus strand or minus strand of the reference sequence.  
		read_seq_aligned is the aligned sequence read, which always in the same orientation of the same read in the corresponding fastq file. 


	SAM format

       		SAM is a standard format for next-gen sequencing read alignments. The details of the format and examples are available at the links below:	
		1) http://samtools.sourceforge.net/SAM1.pdf
	       	2) http://genome.sph.umich.edu/wiki/SAM		

	BED format

		See the format at UCSC http://genome.ucsc.edu/FAQ/FAQformat.html#format1 
		
		NOTE: both ALN and BED format files use 0-based coordinate system while SAM format uses 1-based coordinate system.

ACKNOWLEDGEMENTS
	I would like to thanks all ART users for their feedback and contributions, especially the users listed below.
      	Richard Nielson,  DNASTAR
	Bruno Nevado, CRAG in UAB  
	Lee Katz, US CDC  

