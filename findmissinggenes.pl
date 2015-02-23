#!/usr/bin/perl 

#	input: orthlist, FASTA file, outfile

use warnings;
use strict;

open ORTH, shift @ARGV;
my %totalgenes;
while (<ORTH>) {
	if ($_ =~ /(\S+)\t(\S+)/) {
		$totalgenes{$1}=1;
		$totalgenes{$2}=1;
	}
}

#my @genes = keys(%totalgenes);
#print $#genes;

close ORTH;

open FASTA, shift @ARGV;
open OUT, ">", shift @ARGV;
while (<FASTA>) {
	if ($_ =~ /\>(\S+)/) {
		if (defined ($totalgenes{$1})) {
			next;
		}
		else {
			print OUT $1, "\n";
		}
	}
}

close FASTA;

close OUT;