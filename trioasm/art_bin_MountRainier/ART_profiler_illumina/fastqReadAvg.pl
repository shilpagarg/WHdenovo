#!/usr/bin/perl
# A perl program to retrieve quality scores for each base at a given position
# from a fastq file denoted by the user as a command-line argument.
# @author Jason Myers

use FileHandle;
use strict;
#use warnings;

# Making sure that the correct number of arguments are given by the user and 
# that the program displays an accurate message in case the user does not.

my $numArgs = $#ARGV + 1;

unless( $numArgs == 1) {
    
    print "Usage: perl fastqReadAvg.pl FastqFile\n";
    exit;
}


# Getting the filename to the file of interest
my $infile = $ARGV[0];
my $fh = FileHandle->new;
# Opening the file to be read into and the file to hold the programs output.
#open INFILE, "$infile", or die $!;
if($infile!~/\.gz$/){
       	$fh->open ("<$infile") or die $!;
}	
else{ 
	$fh->open("zcat $infile|") or die $!;
}

my $outFile = $infile . ".txt";
$outFile =~s/^.*\///; 

open OUTFILE, ">>$outFile", or die $!;


# Iterator to determine which line to take from the file.
my $iter = 0;

# variable to switch between the operations associated with the called bases and the quality scores
my $whichLine = 1;

# array to store the current line of called bases
my @BASES = ();

# array to store the current line of quality scores
my @SCORES = ();

# variable to hold initial read of the called bases
my $line1 = '';

# variable to hold initial read of quality scores
my $line2 = '';

# Initializing the main data array that is to be base call at a specific position
# by quality score.
my @data = ();

# Variable holding the number of scores possible (0-39)
my $numScores = 71;

# Variable signifying the number of possible bases (ATGC)
my $numBases = 5;

# Variable holding the largest length of read that is in the file
# initialy 0
my $maxLength = 0;

# The ascII value of '!' in decimal format. Used to normalize quality scores
my $normal = 33;

# Variable associated with the position of 'TGC' base scores in the output file
my $Tpos = 1;
my $Gpos = 2;
my $Cpos = 3;
my $Npos = 4;

# Final data matrix
my @finDat = ();

# Matrix to hold the average scores 
my @avgScore = ();

# Array holding the different quality scores
my @allScores = ();

# Fill the array.
for(my $r = 0; $r < $numScores; $r++){
	$allScores[$r] = $r;
}

# Main while loop that iterates over the fastq file	
#while(<INFILE>) {
while(<$fh>) {
    
    # The lines of interest are the odd lines
    if( ($iter % 2) != 0 ) {
	# If the current line is a called bases line
	if($whichLine){
	    # read in the line
	    $line1 = $_;
	    
	    # remove the new-line character
	    chomp $line1;
	    
	    # split the string up into its individual characters
	    @BASES = split(//, $line1);				

		# record the longest read length to determine how large the data 
		# matrix has grown later
		if($#BASES > $maxLength) {
			# $maxLength = $#BASES; 
			$maxLength = scalar @BASES; # Thanks Thomas Hack for the bug fix  02/09/2014 //weichun
		}
	    
	    # Set the which flag to false because the next odd line will be a score line
	    $whichLine--;
	    
	} else {
	    # read in the line
	    $line2 = $_;
	    
	    # remove the new-line character
	    chomp $line2;
	    
	    # split the line up into individual characters
	    @SCORES = split(//, $line2);
	    
	    # Iterate over the posiions in the score line
	    foreach my $pos1 (0 .. $#SCORES) {
		# Determine which base has been called and it to determine which column of the data
		# matrix needs to be manipulated. The ascII value of the score line associated with
		# this base minus 32 is then used to determine which row of the data matrix needs
		# manipulation. Finally the count is increased for that position in the matrix.
			if( $BASES[$pos1] eq 'A'){
		    	$data[(ord($SCORES[$pos1])-$normal)][($pos1*$numBases)] += 1;			
		    
			} elsif( $BASES[$pos1] eq 'T'){
		   		$data[(ord($SCORES[$pos1])-$normal)][(($pos1*$numBases)+$Tpos)] += 1;			
		    
			} elsif( $BASES[$pos1] eq 'G'){
		    	$data[(ord($SCORES[$pos1])-$normal)][(($pos1*$numBases)+$Gpos)] += 1;			
		    
			} elsif( $BASES[$pos1] eq 'C'){
		   		$data[(ord($SCORES[$pos1])-$normal)][(($pos1*$numBases)+$Cpos)] += 1;			
		    
			} else {
		   		$data[(ord($SCORES[$pos1])-$normal)][(($pos1*$numBases)+$Npos)] += 1;			
			
			}
	    }
	    
	    
	    # the next odd line will be a called bases line so set the which flag back to true
	    $whichLine++;
	}
    }
    
    
    # Increment the iterator
    $iter++;
    
}

# Write the longest length of read to the output file.
print OUTFILE "Longest Length Read: $maxLength\n";

# Put column headings on the file tab delimited
foreach my $pos2 (0 .. ($maxLength - 1)) {
    print OUTFILE "A$pos2\t";

}

foreach my $pos3(0 .. ($maxLength - 1)){
	print OUTFILE "T$pos3\t";

}

foreach my $pos4(0 .. ($maxLength - 1)){
	print OUTFILE "G$pos4\t";

}

foreach my $pos5(0 .. ($maxLength - 1)){
	print OUTFILE "C$pos5\t";

}

foreach my $pos6(0 .. ($maxLength - 1)){
	print OUTFILE "N$pos6\t";
}
		
print OUTFILE "\n";

# Re-arrange the elements of the matrix and calculate the average quality score per position.
for(my $h = 0; $h < $numBases; $h++){
	for(my $t = 0; $t < $numScores; $t++) {
	    for(my $g = 0; $g < ($maxLength); $g++) {
		
			# Set null elements to 0
			if($data[$t][(($g*$numBases)+$h)] < 1){
				$data[$t][(($g*$numBases)+$h)] = 0;
			}			

			# Do the re-arranging 
			$finDat[$t][($g+(($maxLength)*$h))] = $data[$t][(($g*$numBases)+$h)];

			# Do the summations needed to find avg quality scores
			$avgScore[0][($g+(($maxLength)*$h))] += $allScores[$t]*$finDat[$t][($g+(($maxLength)*$h))];
			$avgScore[1][($g+(($maxLength)*$h))] += $finDat[$t][($g+(($maxLength)*$h))];
		
		}
	}
}
	

# Write the data matrix to the output file tab delimited 
for(my $o = 0; $o < ($numScores + 1); $o++) {
    for(my $p = 0; $p < ($maxLength * $numBases); $p++) {

		if($o == $numScores){
			# Make sure the numerators and denominators are not 0.
			if(($avgScore[0][$p] != 0) && ($avgScore[1][$p] != 0)){
				print OUTFILE ($avgScore[0][$p] / $avgScore[1][$p]), "\t"; 
			} else {
				print OUTFILE "0\t";
			}
		} else {
			print OUTFILE "$finDat[$o][$p]\t";
		}

    }
    print OUTFILE "\n";
}	

# Write the average quality score information to the file
for( my $n = 0; $n < 2; $n++){
	for( my $m = 0; $m < ($maxLength *$numBases); $m++){
		print OUTFILE $avgScore[$n][$m], "\t";
	}
	print OUTFILE "\n";
}

# Close the open files.
close(OUTFILE);
close(INFILE);

# Exit the program
exit;
