#!/usr/bin/perl
# A perl program to compute the summation of frequency amongst multiple files for each  
# quality score.
# @author Jason Myers

use strict;

# Print an error message unless the atleast two fastqReadAvg files and an output file are specified
my $numArgs = $#ARGV + 1;

unless( $numArgs >= 3){
    print "Usage: perl summation.pl fastqReadAvg_file1 fastqReadAvg_file2 [fastReadqAvg_fileN] outputFile \n";
    exit;
}


# Arrays to hold the different quality score frequencies by base 
my @COMB_FREQA = ();
my @COMB_FREQT = ();
my @COMB_FREQG = ();
my @COMB_FREQC = ();
my @COMB_FREQN = ();

# A variable to track the maximum length read seen
my $maxLength = 0;

# Loop over the input files
for(my $iter = 0; $iter < ($numArgs - 1); $iter++){

	# Open the file in question
	my $infile = $ARGV[$iter];

	open INFILE, "$infile", or die $!;

	# Set the line counter to 0
	my $count = 0;

	# An array to temporarily store each files frequency lines
	my @TEMP = ();
	
	# A variable to hold the initial grab of the current line
	my $line = '';

	# variable corresponding to the current length
	my $curLength;

	# Loop over each line of the file
	while(<INFILE>){

		if( $count >= 2 && $count <= 73 ){

			# get the line
			$line = $_;
	
			# set the current length to 0.
			$curLength = 0;


			# split the line up into its elements
			@TEMP = split('\t', $line);

			# calculate the current length of read being seen
			$curLength = $#TEMP / 5;
			
			# if the current length is longer than the current maximum length seen
			# set the max length to the current length
			if($curLength > $maxLength){
				$maxLength = $curLength;
				}

			# iterate over the elements of the current line and add them to their corresponging positions
			for(my $iter = 0; $iter < 5; $iter++){
				for(my $k = 0; $k < $curLength; $k++){
					if( $iter == 0){
						$COMB_FREQA[($count - 2)][$k] += $TEMP[(($iter * $curLength) + $k)];

					} elsif( $iter == 1){
						$COMB_FREQT[($count - 2)][$k] += $TEMP[(($iter * $curLength) + $k)];

					} elsif( $iter == 2){
						$COMB_FREQG[($count - 2)][$k] += $TEMP[(($iter * $curLength) + $k)];

					} elsif( $iter == 3){
						$COMB_FREQC[($count - 2)][$k] += $TEMP[(($iter * $curLength) + $k)];

					} else{
						$COMB_FREQN[($count - 2)][$k] += $TEMP[(($iter * $curLength) + $k)];
						
					}
				}
			}

			# set the temp variable back to null
			@TEMP = ();
		}

		# increase the line counter
		$count++;
	}

	# close the current infile
	close(INFILE);
}

# open the outfile
my $outfile = $ARGV[$numArgs - 1];

open OUTFILE, ">>$outfile", or die $!;


# set the first base to 'A'
my $base = 'A';

# write the header out to the output file
for( my $iter = 0; $iter < 5; $iter++){

	if( $iter == 0){
	} elsif( $iter == 1 ){
		$base = 'T';
	} elsif( $iter == 2 ){
		$base = 'G';
	} elsif( $iter == 3 ){
		$base = 'C';
	} else {
		$base = 'N';
	}

	foreach my $pos1 (0 .. ($maxLength - 1)){

		print OUTFILE $base, $pos1, "\t";
	
	}
}

print OUTFILE "\n";

# write the contents of each combined frequency array out to the output
# file in the correct order
for(my $i = 0; $i < 71; $i++){
	for(my $u = 0; $u  < 5; $u++){
		for(my $j = 0; $j < $maxLength; $j++){
			if($u == 0){
				print OUTFILE $COMB_FREQA[$i][$j], "\t"

			} elsif($u == 1){
				print OUTFILE $COMB_FREQT[$i][$j], "\t"

			} elsif($u == 2){
				print OUTFILE $COMB_FREQG[$i][$j], "\t"

			} elsif($u == 3){
				print OUTFILE $COMB_FREQC[$i][$j], "\t"

			} else{
				print OUTFILE $COMB_FREQN[$i][$j], "\t"
			}			
		}
	}
	
	print OUTFILE "\n";
}

# close the output file
close(OUTFILE);

# end  summation.pl
exit;
