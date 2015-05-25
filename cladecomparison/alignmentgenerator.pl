#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;

my (@groups, @reps);

GetOptions(
"g=s@"	=> \@groups,
"r=s@"	=> \@reps);

my @querygenes;

my @comparison = `ls $groups[0]/*.faa`;
foreach my $comp (@comparison) {
	system "esl-sfetch --index $comp";
}

foreach my $num (1..$#groups) {
	my $group = $groups[$num];
	my $rep = $reps[$num];

	open IN, "$reps[0]-$reps[$num].siggenes";
	while (<IN>) {
		chomp;
		my @split = split(/\s+/, $_);
		push @querygenes, $split[0];
	}
	close IN;

	my @members = `ls $group/*.faa`;
	foreach my $member (@members) {
		system "esl-sfetch --index $member";
	}
}

foreach my $geneofinterest (@querygenes) {
	open OUT, "> queryfile.faa";
	foreach my $group (0..$#groups) {
		my @members = `ls $groups[$group]/*.faa`;
		foreach my $file (@members) {
			fetchseq($geneofinterest, $file, $groups[$group], $reps[$group]);
		}
	}
	system "mafft queryfile.faa > alignments/$geneofinterest.afa";
}

sub fetchseq {								# input is fetchseq(gene, faa file, group name)
	my $geneofinterest = $_[0];
	my $file = $_[1];
	my $group = $_[2];
	my $rep = $_[3];
	my $refname;
	if (-e "$reps[0]-$rep/orthlist.dbs") {
		open ORTHLIST, "$reps[0]-$rep/orthlist.dbs";
		while (<ORTHLIST>) {
			if ($_ =~ /(\S+)\s+(\S+)/) {
				if ($1 eq $geneofinterest) {
					$refname = $2;
				}
			}
		}
	}
	if (!$refname) {
		$refname = $geneofinterest;
	}
	print "\n\n$refname\n\n";
	
	if ($file =~ /(.+\/.+).faa/) {
		my $fileID = $1;
		if (-e "$fileID.orths") {
			open ORTHS, "$1.orths";
			while (<ORTHS>) {
				if ($_ =~ /(\S+)\s+(\S+)/) {
					if ($1 eq $refname) {
						print OUT `esl-sfetch -n "$group-$2" $fileID.faa "$2"` or die "failed to find seq $!";
					}
				}
			}
			close ORTHS;
		}
		else {
			print OUT `esl-sfetch -n "$group-$refname" $fileID.faa "$refname"`;
		}
	}
}
	

