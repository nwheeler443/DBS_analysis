#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# could use input of folder names to name groups?

my @inputfiles = `ls */*/*.scan.t3ees`;

# get filtered list for each isolate from each group

my %scores;

foreach my $input (@inputfiles) {
	my $category;
	my $isolate;
	if ($input =~ /(.+)\/.+\/(.+).scan.t3ees/) {
		$category = $1;
		$isolate = $2;
	}
	open IN, $input;
	while (<IN>) {
		if ($_ =~ /#/) {
			next;
		}
		else {
			my @split = split(/\s+/, $_);
			if (defined($scores{$category}{$isolate}{$split[3]}{model})) {
				if ($scores{$category}{$isolate}{$split[3]}{score} < $split[7]) {
					$scores{$category}{$isolate}{$split[3]}{model} = $split[0];
					$scores{$category}{$isolate}{$split[3]}{score} = $split[7];
				}
			}
			else {
				$scores{$category}{$isolate}{$split[3]}{model} = $split[0];
				$scores{$category}{$isolate}{$split[3]}{score} = $split[7];
			}
		}
	}
}
close IN;

my @effectors;

open OUT, "> T3SS_Genes.txt";
foreach my $category (keys(%scores)) {
	foreach my $species (keys(%{$scores{$category}})) {
		foreach my $gene (keys(%{$scores{$category}{$species}})) {
			print OUT "$category\t";
			print OUT "$species\t";
			print OUT "$gene\t";
			print OUT "$scores{$category}{$species}{$gene}{'model'}\t";
			print OUT "$scores{$category}{$species}{$gene}{'score'}\n";
			if ($scores{$category}{$species}{$gene}{'model'} ~~ @effectors) {
				next;
			}
			else {
				push @effectors, $scores{$category}{$species}{$gene}{'model'};
			}
		}
	}
}

my $found = 0;
foreach my $category (keys(%scores)) {
	open OUT, "> T3SS_scores_$category.txt";
	foreach my $effector (@effectors) {
		print OUT "$effector";
		foreach my $species (keys(%{$scores{$category}})) {
			my $bestscore = 0;
			foreach my $gene (keys(%{$scores{$category}{$species}})) {
				if ($scores{$category}{$species}{$gene}{model} eq $effector) {
					$found = 1;
					if ($scores{$category}{$species}{$gene}{score} > $bestscore) {
						$bestscore = $scores{$category}{$species}{$gene}{score};
					}
				}
				
			}
			if ($found ==1) {
				print OUT "\t$bestscore";
			}
			if ($found ==0) {
				print OUT "\tNA";
			}
			$found = 0;
		}
		print OUT "\n";
	}
}
