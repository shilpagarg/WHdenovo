open(F1, $ARGV[0]) or die $!;
open(F2, $ARGV[1]) or die $!;
while(<F1>){
     	chomp;
       	@aa=split /\s+/;
       	$aa[0]=($aa[0]-1)*2;
       	printf "%s\n", join "\t", @aa;
	my $bf=<F2>;
     	chomp;
       	@aa=split /\s+/,$bf; $aa[0]=($aa[0]-1)*2+1; printf "%s\n", join "\t", @aa;
}
