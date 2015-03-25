#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;

my $seqfile = shift @ARGV;

open IN, $seqfile;
my @sequence;

while (<IN>) {
	chomp;
	if ($_ =~ />/) {
		next;
	}
	else {
		push @sequence, split("", $_);
	}
}

close IN;

open OUT, "> $seqfile.numbered";
foreach my $res (0..$#sequence) {
	print OUT $res+1, "\t$sequence[$res]\n";
}

close OUT;