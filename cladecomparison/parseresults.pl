#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my @keys = ("Pathogenic", "Rhizosphere", "Environmental");
my @refs = ("Pathogenic/Pto_DC3000_P/Pto_DC3000_P_AE16853", "Rhizosphere/Pfl_PCL1751_R/Pfl_PCL1751_R_CP010896", "Environmental/Pfl_Pf0-1_E/Pfl_Pf0-1_E_NC_007492");
my @dbsfiles = ("path-rhiz.dbs/results.dbs", "path-rhiz.dbs/results.dbs", "path-env.dbs/results.dbs");		# will only get pathogen genes with an ortholog in rhizosphere dbs comparison

foreach my $comparison (0..2) {
	my $key = $keys[$comparison];
	my $ref = $refs[$comparison];
	print "$ref\n";
	my @scanfiles = `ls $key/*/*.scan`;
	my %orths;
	my $archfile = $dbsfiles[$comparison];
	my %genearchs;
	# makes a hash with gene architectures for each gene and makes orthlist for the reference species
	my $index;
	if ($comparison == 0) {$index = 0}
	else {$index = 1}
	open DBSFILE, "$archfile";
	while (<DBSFILE>) {
		if ($_ =~/#/) {next;}
		my @split = split(/\s+/, $_);
		push @{$genearchs{$split[$index]}{domains}}, $split[2];
		push @{$genearchs{$split[$index]}{scores}}, $split[9];
		$orths{$ref}{$split[$index]} = $split[$index];
	}

	# fills out the rest of the ortholog hash
	foreach my $scan (@scanfiles) {
		if ($scan =~ /($key\/.+\/.+).scan/) {
			open ORTHS, "$1.orths";
			my $name = $1;
			while (<ORTHS>) {
				if ($_ =~ /(\S+)\s+(\S+)/) {
					$orths{$name}{$2} = $1;
				}
			}
		}
	}
	
	my %fullscores;
	my %scores;
	my @names;
	
	# makes a hash of scores for each gene
	foreach my $scan (@scanfiles) {
		my $name;
		if ($scan =~ /($key\/.+\/.+).scan/) {
			$name = $1;
		}
		push @names, $name;
		open SCAN, $scan;
		while (<SCAN>) {
			if ($_ =~ /#/) {next;}
			my @split = split(/\s+/, $_);
			if (defined($orths{$name}{$split[3]})) {		# if there is a representative ortholog for this gene
					push @{$fullscores{$name}{$split[3]}{$split[1]}}, $split[7];
			}
		}
		close SCAN;
		# now need to go through scores hash and only keep the ones that match the domain architecture (remove multiples)
		foreach my $gene (keys(%{$fullscores{$name}})) {
			foreach my $domain (@{$genearchs{$orths{$name}{$gene}}{domains}}) {
				push @{$scores{$name}{$gene}}, $fullscores{$name}{$gene}{$domain}[0];
				shift @{$fullscores{$name}{$gene}{$domain}};
			}
		}
	}
	
	my %scoresums;
	# initialises score list with NAs, will replace them if an ortholog is present
	foreach my $refgene (keys(%{$orths{$ref}})) {
		foreach my $pos (0..$#names) {
			push @{$scoresums{$refgene}}, "NA";
		}
	}
	# subs in ortholog scores
	foreach my $pos (0..$#names) {
		my $name = $names[$pos];				# for each member of the group
		foreach my $gene (keys(%{$scores{$name}})) {		# for each gene
			my $runningtotal = 0;
				foreach my $score (@{$scores{$name}{$gene}}) {			# why won't this work?! still getting a valid output file
					$runningtotal += $score;							# add up the individual domain scores
				}
			$scoresums{$orths{$name}{$gene}}[$pos] = $runningtotal;
			print $runningtotal;
		}
	}
	
	open OUT, "> $key.scores.NAs.txt";
	foreach my $pos (0..$#names) {
		if ($names[$pos] =~ /^$key\/(.+)\/.+$/) {
			print OUT "\t$1";
		}
	}
	print OUT "\n";
	foreach my $refgene (keys(%scoresums)) {
		print OUT $refgene;
		foreach my $pos (0..$#names) {
			print OUT "\t", $scoresums{$refgene}[$pos];
		}
		print OUT "\n";
	}
	close OUT;
}
