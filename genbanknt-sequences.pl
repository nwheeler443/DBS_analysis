#!/usr/bin/perl

## turns genbank entries into nucleotide sequences

use warnings;
use strict;
use Data::Dumper;

my $genbank = shift @ARGV;
my $fasta = shift @ARGV;

system "esl-sfetch --index $fasta";

open GB, $genbank;
open FASTA, $fasta;

my $genomename;
while (<FASTA>) {
	if ($_ =~ />(.*)/) {
		$genomename = $1;
	}
}

my $coordinates;
my $genename;
my $printme1 = 0;
my $printme2 = 0;

my $outfile = shift @ARGV;
open OUT, ">> $outfile";

# add a line to delete the old .faa file

while (<GB>) {
	if ($_ =~ /CDS\s+(\d+\.\.\d+)/) {
		$coordinates = $1;
		$printme1 = 1;
	}
	if ($_ =~ /CDS\s+complement\((\d+\.\.\d+)\)/) {
		$coordinates = $1;
		$printme1 = 1;
	}
	if ($_ =~ /\/locus_tag="(\S+)"/) {
		$genename = $1;
		$printme2 = 1;
	}
	if ($printme1 == 1 & $printme2 == 1) {
		my $sequence = `esl-sfetch -c $coordinates -n $genename $fasta $genomename`;
		print OUT "$sequence";
		$printme1 = 0;
		$printme2 = 0;
	}
}