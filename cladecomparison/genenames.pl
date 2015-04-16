#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my @genesofinterest;
open PATHENV, "path-env.NAs.siggenes";
while (<PATHENV>) {
	chomp;
	push @genesofinterest, $_;
}

open FAA, "Pathogenic/Pto_DC3000_P/Pto_DC3000_P_AE16853.faa";
open OUT, "> path-env.NAs.genes";
while (<FAA>) {
	if ($_ =~ />(\S+)\s+.+protein=(.+)\]\s+\[protein/) {
		if ($1 ~~ @genesofinterest) {
			print OUT "$1\t$2\n";
		}
	}
}

my @genesofinterest2;
open PATHENV, "path-rhiz.NAs.siggenes";
while (<PATHENV>) {
	chomp;
	push @genesofinterest2, $_;
}

open FAA, "Pathogenic/Pto_DC3000_P/Pto_DC3000_P_AE16853.faa";
open OUT, "> path-rhiz.NAs.genes";
while (<FAA>) {
	if ($_ =~ />(\S+)\s+.+protein=(.+)\]\s+\[protein/) {
		if ($1 ~~ @genesofinterest2) {
			print OUT "$1\t$2\n";
		}
	}
}

my @genesofinterest;
open PATHENV, "path-env.siggenes";
while (<PATHENV>) {
	chomp;
	push @genesofinterest, $_;
}

open FAA, "Pathogenic/Pto_DC3000_P/Pto_DC3000_P_AE16853.faa";
open OUT, "> path-env.genes";
while (<FAA>) {
	if ($_ =~ />(\S+)\s+.+protein=(.+)\]\s+\[protein/) {
		if ($1 ~~ @genesofinterest) {
			print OUT "$1\t$2\n";
		}
	}
}

my @genesofinterest2;
open PATHENV, "path-rhiz.siggenes";
while (<PATHENV>) {
	chomp;
	push @genesofinterest2, $_;
}

open FAA, "Pathogenic/Pto_DC3000_P/Pto_DC3000_P_AE16853.faa";
open OUT, "> path-rhiz.genes";
while (<FAA>) {
	if ($_ =~ />(\S+)\s+.+protein=(.+)\]\s+\[protein/) {
		if ($1 ~~ @genesofinterest2) {
			print OUT "$1\t$2\n";
		}
	}
}
