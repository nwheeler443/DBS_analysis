#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my @pathogenic = `ls Pathogenic/*/*.scan`;
my @environmental = `ls Environmental/*/*.scan`;
my @rhizosphere = `ls Rhizosphere/*/*.scan`;

my @querygenes;

open IN, "path-env.NAs.siggenes";
while (<IN>) {
	chomp;
	push @querygenes, $_;
}

my $comparison = "env";

foreach my $geneofinterest (@querygenes) {
	open OUT, "> queryfile.faa";
	
	foreach my $path (@pathogenic) {
		if ($path =~ /(Pathogenic\/.+\/.+).scan/) {
			my $fileID = $1;
			#system "esl-sfetch --index $fileID.faa";
			if (-e "$fileID.orths") {
				open ORTHS, "$1.orths";
				while (<ORTHS>) {
					if ($_ =~ /(\S+)\s+(\S+)/) {
						if ($1 eq $geneofinterest) {
							print OUT `esl-sfetch $fileID.faa "$2"` or die "failed to find seq $!";
						}
					}
				}
				close ORTHS;
			}
			else {
				print OUT `esl-sfetch $fileID.faa "$geneofinterest"`;
			}
		}
	}
	
	print OUT "#####\n";
	
	if ($comparison eq "env") {
		my $refname;
		open PATHENV, "path-env.dbs/orthlist.dbs";
		while (<PATHENV>) {
			if ($_ =~ /(\S+)\s+(\S+)/) {
				if ($1 eq $geneofinterest) {
					$refname = $2;
				}
			}
		}
		foreach my $env (@environmental) {
			if ($env =~ /(Environmental\/.+\/.+).scan/) {
				my $fileID = $1;
				system "esl-sfetch --index $fileID.faa";
				if (-e "$fileID.orths") {
					open ORTHS, "$1.orths";
					while (<ORTHS>) {
						if ($_ =~ /(\S+)\s+(\S+)/) {
							if ($1 eq $refname) {
								print OUT `esl-sfetch $fileID.faa "$2"` or die "failed to find seq $!";
							}
						}
					}
					close ORTHS;
				}
				else {
					print OUT `esl-sfetch $fileID.faa "$refname"`;
				}
			}
		}
		
	}
	elsif ($comparison eq "rhiz") {
		my $refname;
		open PATHRHIZ, "path-rhiz.dbs/orthlist.dbs";
		while (<PATHRHIZ>) {
			if ($_ =~ /(\S+)\s+(\S+)/) {
				if ($1 eq $geneofinterest) {
					$refname = $2;
				}
			}
		}
		foreach my $rhiz (@rhizosphere) {
			if ($rhiz =~ /(Rhizosphere\/.+\/.+).scan/) {
				my $fileID = $1;
				system "esl-sfetch --index $fileID.faa";
				if (-e "$fileID.orths") {
					open ORTHS, "$1.orths" or die "couldn't open file $1.orths $!";
					my $fileID = $1;
					while (<ORTHS>) {
						if ($_ =~ /(\S+)\s+(\S+)/) {
							if ($1 eq $refname) {
								print OUT `esl-sfetch $fileID.faa "$2"` or die "failed to find seq $!";
							}
						}
					}
					close ORTHS;
				}
				else {
					print OUT `esl-sfetch $fileID.faa "$refname"`;
				}
			}
		}
	}
	
	my $filename = $geneofinterest;
	if ($geneofinterest =~ /lcl\|(\S+)/) {
		$filename = $1;
	}
	print "\n\n$filename\n\n";
	system "mafft queryfile.faa > alignments/$comparison-$filename.afa";
}


@querygenes = ();

$comparison = "rhiz";
open IN, "path-rhiz.NAs.siggenes";
while (<IN>) {
	chomp;
	push @querygenes, $_;
}

foreach my $geneofinterest (@querygenes) {
	open OUT, "> queryfile.faa";
	
	foreach my $path (@pathogenic) {
		if ($path =~ /(Pathogenic\/.+\/.+).scan/) {
			my $fileID = $1;
			#system "esl-sfetch --index $fileID.faa";
			if (-e "$fileID.orths") {
				open ORTHS, "$1.orths";
				while (<ORTHS>) {
					if ($_ =~ /(\S+)\s+(\S+)/) {
						if ($1 eq $geneofinterest) {
							print OUT `esl-sfetch $fileID.faa "$2"` or die "failed to find seq $!";
						}
					}
				}
				close ORTHS;
			}
			else {
				print OUT `esl-sfetch $fileID.faa "$geneofinterest"`;
			}
		}
	}
	
	print OUT "#####\n";
	
	if ($comparison eq "env") {
		my $refname;
		open PATHENV, "path-env.dbs/orthlist.dbs";
		while (<PATHENV>) {
			if ($_ =~ /(\S+)\s+(\S+)/) {
				if ($1 eq $geneofinterest) {
					$refname = $2;
				}
			}
		}
		foreach my $env (@environmental) {
			if ($env =~ /(Environmental\/.+\/.+).scan/) {
				my $fileID = $1;
				system "esl-sfetch --index $fileID.faa";
				if (-e "$fileID.orths") {
					open ORTHS, "$1.orths";
					while (<ORTHS>) {
						if ($_ =~ /(\S+)\s+(\S+)/) {
							if ($1 eq $refname) {
								print OUT `esl-sfetch $fileID.faa "$2"` or die "failed to find seq $!";
							}
						}
					}
					close ORTHS;
				}
				else {
					print OUT `esl-sfetch $fileID.faa "$refname"`;
				}
			}
		}
		
	}
	elsif ($comparison eq "rhiz") {
		my $refname;
		open PATHRHIZ, "path-rhiz.dbs/orthlist.dbs";
		while (<PATHRHIZ>) {
			if ($_ =~ /(\S+)\s+(\S+)/) {
				if ($1 eq $geneofinterest) {
					$refname = $2;
				}
			}
		}
		foreach my $rhiz (@rhizosphere) {
			if ($rhiz =~ /(Rhizosphere\/.+\/.+).scan/) {
				my $fileID = $1;
				system "esl-sfetch --index $fileID.faa";
				if (-e "$fileID.orths") {
					open ORTHS, "$1.orths" or die "couldn't open file $1.orths $!";
					my $fileID = $1;
					while (<ORTHS>) {
						if ($_ =~ /(\S+)\s+(\S+)/) {
							if ($1 eq $refname) {
								print OUT `esl-sfetch $fileID.faa "$2"` or die "failed to find seq $!";
							}
						}
					}
					close ORTHS;
				}
				else {
					print OUT `esl-sfetch $fileID.faa "$refname"`;
				}
			}
		}
	}
	
	my $filename = $geneofinterest;
	if ($geneofinterest =~ /lcl\|(\S+)/) {
		$filename = $1;
	}
	print "\n\n$filename\n\n";
	system "mafft queryfile.faa > alignments/$comparison-$filename.afa";
}
