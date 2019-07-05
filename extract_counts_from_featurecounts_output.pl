#!/usr/bin/perl -w
use strict;
my @file = glob "*.tsv";
print "@file\n";
my $dir = `pwd`;
system("mkdir results");
print "$dir\n";
foreach my $a (@file) {
	open IN, "$a" or die "$!";
	open OUT,">$a.txt" or die "$!";

	while (my $line = <IN>) {
		chomp $line;
	        next if ($line =~ /^#/);
		my ($geneid,$count) = (split /\t/,$line)[0,6];
			
		print OUT "$geneid\t$count\n";
	
	}
	close OUT;
	close IN;
}
