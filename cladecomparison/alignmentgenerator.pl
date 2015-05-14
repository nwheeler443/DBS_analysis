#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my @pathogenic = `ls Pathogenic/*/*.scan`;
my @environmental = `ls Environmental/*/*.scan`;
my @rhizosphere = `ls Rhizosphere/*/*.scan`;

foreach my $file (@pathogenic) {
	eslindex($file);
}
foreach my $file (@environmental) {
	eslindex($file);
}
foreach my $file (@rhizosphere) {
	eslindex($file);
}
sub eslindex {
	my $file = $_[0];
	if ($file =~ /(.+\/.+\/.+).scan/) {
		my $fileID = $1;
		system "esl-sfetch --index $fileID.faa";
	}
}

sub aligngenes {
	my $geneofinterest = $_[0];
	open OUT, "> queryfile.faa";
	
	foreach my $path (@pathogenic) {
		if ($path =~ /(Pathogenic\/.+\/.+).scan/) {
			my $fileID = $1;
			if (-e "$fileID.orths") {
				open ORTHS, "$1.orths";
				while (<ORTHS>) {
					if ($_ =~ /(\S+)\s+(\S+)/) {
						if ($1 eq $geneofinterest) {
							print OUT `esl-sfetch -n "path-$1" $fileID.faa "$2"` or die "failed to find seq $!";
						}
					}
				}
				close ORTHS;
			}
			else {
				print OUT `esl-sfetch -n "path-$1" $fileID.faa "$geneofinterest"`;
			}
		}
	}
	
	print OUT "#####\n";
	
	foreach my $file (@rhizosphere) {
		fetchseq($geneofinterest, $file, "rhiz");
	}
	foreach my $file (@environmental) {
		fetchseq($geneofinterest, $file, "env");
	}
	
	sub fetchseq {
		my $geneofinterest = $_[0];
		my $file = $_[1];
		my $prefix = $_[2];
		my $refname;
		open ORTHLIST, "path-$prefix.dbs/orthlist.dbs";
		while (<ORTHLIST>) {
			if ($_ =~ /(\S+)\s+(\S+)/) {
				if ($1 eq $geneofinterest) {
					$refname = $2;
				}
			}
		}
		if ($file =~ /(.+\/.+\/.+).scan/) {
			my $fileID = $1;
			if (-e "$fileID.orths") {
				open ORTHS, "$1.orths";
				while (<ORTHS>) {
					if ($_ =~ /(\S+)\s+(\S+)/) {
						if ($1 eq $refname) {
							print OUT `esl-sfetch -n "$prefix-$1" $fileID.faa "$2"` or die "failed to find seq $!";
						}
					}
				}
				close ORTHS;
			}
			else {
				print OUT `esl-sfetch -n "$prefix-$refname" $fileID.faa "$refname"`;
			}
		}
	}
	
	my $filename;
	if ($geneofinterest =~ /lcl\|(\S+)/) {
		$filename = $1;
	}
	print "\n\n$filename\n\n";
	system "mafft queryfile.faa > alignments/$filename.afa";
}

my @querygenes;
open IN, "path-env.NAs.siggenes";
while (<IN>) {
	chomp;
	push @querygenes, $_;
}
close IN;
open IN, "path-rhiz.NAs.siggenes";
while (<IN>) {
	chomp;
	if ($_ ~~ @querygenes) {
		next;
	}
	else {
		push @querygenes, $_;
	}
}
close IN;

foreach my $geneofinterest (@querygenes) {
	open OUT, "> queryfile.faa";
	aligngenes($geneofinterest);
}


