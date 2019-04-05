#!/usr/bin/perl
#convert ART's MAP format to the UCSC BED format
#weichun huang at whduke@gmail.com
use strict;
my @mapFile=@ARGV;
if(@ARGV<2){
	print "map2bed converts ART's map files to a BED file\n\n";
	print "USAGE: ";
	print "$0 out_bed_file.bed in_map_file_1 [ in_map_file_2 ...]\n\n";
	exit(1);
}

my $score=500; #all scores are 500.
#bed file: ref_id st_pos end_pos read_id 500 strand
my $bedFile=shift @mapFile;
open(BED, ">$bedFile") or die $!;
my %reflen;
my $readlen=0;
foreach my $af (@mapFile){
	open(MAP, "<$af") or die $!;
	my $i=0;
	#read header
	while (<MAP>){
		chomp;
		if($i==0){
			if(!/^##ART_SOLiD\t\S+\t(\d+)/){
				print STDERR "File $af is not an ART_SOLiD MAP file.\n Exit ...\n";
				exit(1);
			}
			$readlen=$1;
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
	while (<MAP>){
		chomp;
		next if(!/^\w/);
	       	my ($ref_id, $read_id, $st_pos, $strand, @aa)=split /\t/;
	       	if($strand eq "+"){
		       	printf BED "%s\t%d\t%d\t%s\t%d\t%s\n", ($ref_id,  $st_pos, $st_pos+$readlen, $read_id, $score, $strand); 
		}
	       	else{
		       	my $stop=$reflen{$ref_id}-$st_pos;
		       	printf BED "%s\t%d\t%d\t%s\t%d\t%s\n", ($ref_id,  $stop-$readlen, $stop, $read_id, $score, $strand); 
		}
       	} 
}

