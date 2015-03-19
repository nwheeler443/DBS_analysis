#!/usr/bin/perl

#want to be able to give a profile HMM and a protein sequence, then create an output file indicating which positions in the protein would best tolerate insertions, deletions and insert extensions.

use warnings;
use strict;
use Data::Dumper;

my $hmmfile = shift @ARGV;
my $alignmentfile = shift @ARGV;

my %embitscores;
my %insprobabilities;
my %transprobabilities;
my $position = 0;
my $key;

open HMM, $hmmfile;
while (<HMM>) {
	if ($_ =~ /^\s+(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+)$/) {
		$key = "pos$position";
		@{$transprobabilities{$key}} = split("  ", $1)
	}
}
close HMM;

# read in alignment file, make an array called @states which has keys for insertions and matches, and put complete sequence into an array

open ALI, $alignmentfile;
my $seqnum = 1;
my $name;
my @firstseq;
my @states;
while (<ALI>) {
	chomp;
	if ($_ =~ /#=GC RF\s+(.+)$/) {
		push @states, split("", $1);
	}
	next if ($_ =~ /#/);
	chomp;
	if ($_ =~ /(\S+)\s+(\S+)$/) {
		push @firstseq, split("", $2);
		$name = $1;
	}
}
close ALI;

