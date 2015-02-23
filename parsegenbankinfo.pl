#!/usr/bin/perl 

use warnings;
use strict;
use Data::Dumper;

my $genbank = shift @ARGV;

open GB, $genbank;

my $genename="";
my $locustag="";
my $note = "";
my $product="";
my $printme = 0;

my $outfile = shift @ARGV;
open OUT, "> $outfile";

# add a line to delete the old .faa file

while (<GB>) {
	if ($_ =~ /^\s+CDS/) {
		$locustag="";
		$product="";
		$genename="";
		$note ="";
	}
	if ($_ =~ /\/gene="(\S+)"/) {
		$genename = $1;
	}
	if ($_ =~ /\/locus_tag="(\S+)"/) {
		$locustag = $1;
	}
	if ($_ =~ /\/note="(.*)"/) {
		$note = $1;
	}
	if ($_ =~ /\/product="(.*)"/) {
		$product = $1;
		$printme = 1;
	}
	if ($printme == 1) {
		print OUT "$locustag\t$genename\t$product\t$note\n";
		$printme = 0;
		
	}
}

close GB;