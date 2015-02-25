#!/usr/bin/perl 

#	input: orthlist, FASTA1, FASTA2, outfile

use warnings;
use strict;
use Data::Dumper;

open ORTH, shift @ARGV;
my %orthlist;
while (<ORTH>) {
	if ($_ =~ /(\S+)\t(\S+)/) {
		$orthlist{$1}=$2;
	}
}

close ORTH;

open FASTA1, shift @ARGV;
my %fasta1;
my $seqname;
while (<FASTA1>) {
	chomp;
	if ($_ =~ />(\w+)/) {
		$seqname=$1;
	}
	else {
		push @{$fasta1{$seqname}}, split("", $_);
	}
}

close FASTA1;

open FASTA2, shift @ARGV;
my %fasta2;
while (<FASTA2>) {
	chomp;
	if ($_ =~ />(\w+)/) {
		$seqname=$1;
	}
	else {
		push @{$fasta2{$seqname}}, split("", $_);
	}
}

close FASTA2;

open OUT, ">", shift @ARGV;
foreach my $seq (keys(%orthlist)) {
	my $otherseq = $orthlist{$seq};
	my $sameseq = 1;
	foreach my $pos (0..$#{$fasta1{$seq}}) {
		if ($fasta1{$seq}[$pos] eq $fasta2{$otherseq}[$pos]) {
			next;
		}
		else {$sameseq = 0}
	}
	if ($sameseq==1) {
		print OUT "$seq\n";
	}
}

close OUT;