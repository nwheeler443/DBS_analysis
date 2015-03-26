#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# could use input of folder names to name groups?

my @pathogenic = `ls Pathogenic/*/*.scan`;
my @environmental = `ls Environmental/*/*.scan`;
my @rhizosphere = `ls Rhizosphere/*/*.scan`;

my $ref = "Environmental/Pfl_Pf0-1_E/Pfl_Pf0-1_E_NC_007492";
#my $ref = "Rhizosphere/Pfl_PCL1751_R/Pfl_PCL1751_R_CP010896";
my $class = "Environmental";

### CHOOSING TO GO WITH A SUMMED SCORE ACROSS DOMAINS FOR SIMPLICITY'S SAKE

my %genearchs;
my %orths;

# domain architectures
open PATHREP, "path-env.dbs/results.dbs";
while (<PATHREP>) {
	if ($_ =~/#/) {
		next;
	}
my @split = split(/\s+/, $_);
	push @{$genearchs{$split[1]}{domains}}, $split[2];
	push @{$genearchs{$split[1]}{scores}}, $split[9];
$orths{$ref}{$split[1]} = $split[1];		# losing a lot of data at this stage since there are only 1700 unique genes consistent across the pethogen-environmental comparison with Pfam domain hits
}

# ortholog lookup table
foreach my $path (@environmental) {
	if ($path =~ /($class\/.+\/.+).scan/) {
		open ORTHS, "$1.orths";
		my $name = $1;
		while (<ORTHS>) {
			if ($_ =~ /(\S+)\s+(\S+)/) {
				$orths{$name}{$2} = $1;
			}
		}
		my @orthnumber = keys(%{$orths{$name}});
	}
}

# score hashes
my %scores;
my @pathnames;

foreach my $path (@environmental) {
	my $name;
	if ($path =~ /($class\/.+\/.+).scan/) {
		$name = $1;
	}
	push @pathnames, $name;
	open SCAN, $path;
	while (<SCAN>) {
		if ($_ =~ /#/) {
			next;
	}
		my @split = split(/\s+/, $_);
	if (defined($orths{$name}{$split[3]})) {
		if ($split[1] ~~ @{$genearchs{$orths{$name}{$split[3]}}{domains}}) {
			push @{$scores{$name}{$split[3]}}, $split[7];
			}
		}
	}
}

# hash of sum of scores for each gene - some of these aren't working
my %scoresums;
foreach my $refgene (keys(%{$orths{$ref}})) {
	foreach my $pos (0..$#pathnames) {
		push @{$scoresums{$refgene}}, 0;
	}
}

foreach my $key (keys(%{$orths{$ref}})) {
	foreach my $pos (0..$#pathnames) {
		my $name = $pathnames[$pos];
		foreach my $gene (keys(%{$scores{$name}})) {
			my $runningtotal = 0;
			foreach my $score (@{$scores{$name}{$gene}}) {
				$runningtotal += $score;
			}
				$scoresums{$orths{$name}{$gene}}[$pos] = $runningtotal;
		}
	}
}

open OUT, "> environmentalscores.txt";
foreach my $refgene (keys(%scoresums)) {
	print OUT $refgene;
	foreach my $pos (0..$#pathnames) {
		print OUT "\t", $scoresums{$refgene}[$pos];
	}
	print OUT "\n";
}
close OUT;