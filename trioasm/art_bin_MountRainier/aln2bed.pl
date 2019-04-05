#!/usr/bin/perl
#convert Aln format to the UCSC BED format
#weichun huang at whduke@gmail.com
use strict;
my @alnFile=@ARGV;
if(@ARGV<2){
	print "aln2bed converts ALN files to a BED file\n\n";
	print "USAGE: ";
	print "aln2bed out_bed_file.bed in_aln_file_1 [ in_aln_file_2 ...]\n\n";
	exit(1);
}

my $score=500; #all scores are 500.
#bed file: ref_id st_pos end_pos read_id 500 strand
my $bedFile=shift @alnFile;
open(BED, ">$bedFile") or die $!;
my %reflen;
foreach my $af (@alnFile){
	open(ALN, "<$af") or die $!;
	my $i=0;
	#read header
	while (<ALN>){
		chomp;
		if($i==0){
			if(!/^##ART/){
				print STDERR "File $af is not an ART ALN file.\n Exit ...\n";
				exit(1);
			}
			$i=1;
		}
		elsif(/^\@SQ\t(\S+)\t(\d+)/){
			$reflen{"$1"}=$2;
		}
		elsif(/^\@CM/){
			next;
		}
		elsif(/^##Header\s+End/){
			last;
		}
	}

	#read aln records 
	while (<ALN>){
		chomp;
		if(s/^>//){
		       	my ($ref_id, $read_id, $st_pos, $strand)=split /\t/;
		       	my $len=0;
		       	while (<ALN>){
				chomp;
			       	if(/^\w/){
				       	s/-//g;
				       	$len=length; #lenth of reference part
			       	}
			       	last;
		       	}
			if($strand eq "+"){
			       	printf BED "%s\t%d\t%d\t%s\t%d\t%s\n", ($ref_id,  $st_pos, $st_pos+$len, $read_id, $score, $strand); 
			}
			else{
				my $stop=$reflen{$ref_id}-$st_pos;	
			       	printf BED "%s\t%d\t%d\t%s\t%d\t%s\n", ($ref_id,  $stop-$len, $stop, $read_id, $score, $strand); 
			}
		}
       	} 
}

