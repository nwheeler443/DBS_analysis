#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;

my (@groups, @reps, @dbsfiles);

GetOptions(
				"g=s@"	=> \@groups,
				"r=s@"	=> \@reps,
				"d=s@"	=> \@dbsfiles);

foreach my $num (0..$#groups) {
	my $group = $groups[$num];
	my $rep = $reps[$num];
	my $repcode = "$group/$rep";
	my @scanfiles = `ls $group/*.scan`;
	my %orths;
	my $archfile = $dbsfiles[$num];
	my %genearchs;
	# makes a hash with gene architectures for each gene and makes orthlist for the reference species
	my $index;
	if ($num == 0) {$index = 0}
	else {$index = 1}
	open DBSFILE, "$archfile";
	while (<DBSFILE>) {
		if ($_ =~/#/) {next;}
		my @split = split(/\s+/, $_);
		push @{$genearchs{$split[$index]}{domains}}, $split[2];
		push @{$genearchs{$split[$index]}{scores}}, $split[9];
		$orths{$repcode}{$split[$index]} = $split[$index];
	}

	# fills out the rest of the ortholog hash
	foreach my $scan (@scanfiles) {
		if ($scan =~ /($group\/.+).scan/) {			# might want to put in a next if $rep type thing here
			my $name = $1;
			if ($name eq $repcode) {
				next;
			}
			open ORTHS, "$1.orths";
			print "$1.orths\n";
			while (<ORTHS>) {
				if ($_ =~ /(\S+)\s+(\S+)/) {
					$orths{$name}{$2} = $1;
				}
			}
			my @genes = keys(%{$orths{$name}});
		}
	}
	close ORTHS;
	
	my %fullscores;
	my %scores;
	my @names;
	
	# makes a hash of scores for each gene
	foreach my $scan (@scanfiles) {
		my $name;
		if ($scan =~ /($group\/.+).scan/) {			# strange that the group is included here...
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
		my @genes = keys(%{$scores{$name}});
		print "\n\n", $#genes, "\n\n";
	}
	
	my %scoresums;
	# initialises score list with NAs, will replace them if an ortholog is present
	foreach my $refgene (keys(%{$orths{$repcode}})) {
		foreach my $pos (0..$#names) {
			push @{$scoresums{$refgene}}, "NA";
		}
	}
	# subs in ortholog scores
	foreach my $pos (0..$#names) {
		my $name = $names[$pos];				# for each member of the group
		foreach my $gene (keys(%{$scores{$name}})) {		# for each gene
			my $runningtotal = 0;
			foreach my $score (@{$scores{$name}{$gene}}) {
				if (defined($score)) {
					$runningtotal += $score or die Dumper (\@{$scores{$name}{$gene}});							# add up the individual domain scores
				}
			}
			$scoresums{$orths{$name}{$gene}}[$pos] = $runningtotal;
		}
	}
	
	open OUT, "> $rep.scores.txt";
	foreach my $pos (0..$#names) {
		if ($names[$pos] =~ /^$group\/(.+)$/) {
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
