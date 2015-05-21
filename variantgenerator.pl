#!/usr/bin/perl

# given the humsavar database and a FASTA file of all the original sequences in the database, will genrate variant sequences for all of them

use warnings;
use strict;
use Data::Dumper;

my %variants;
open IN, "humsavar_variants.txt";			# only includes genes that were scored by provean
while (<IN>) {
	if ($_ =~ /^(\S+)\s+(\S+)/) {
		push @{$variants{$1}}, $2;
	}
}

my @queryproteins = keys(%variants);

open OUT, "> missing_query_sequences.txt";

my @sequences = `ls missed_fasta/*.fasta`;
foreach my $file (@sequences) {
	if ($file =~ /missed_fasta\/(.+).fasta/) {
		if ($1 ~~ @queryproteins) {
			my $protein_id = $1;
			my @sequence;
			open FASTA, "missed_fasta/$1.fasta";			# fasta file
			while (<FASTA>) {
				chomp;
				next if($_ =~ /^\>/);
				push (@sequence, split("", $_));
			}
			close FASTA;
			
			my %variantseqs;
			print OUT ">$1_wild_type\n", join("", @sequence), "\n";
			
			foreach my $variant (@{$variants{$1}}) {
				my $pos;
				my $start;
				my $new;
				if ($variant =~ /^(\w)(\d+)(\w)$/) {
					$start = $1;
					$pos = $2;
					$new = $3;
				}
				@{$variantseqs{$variant}} = @sequence;
				if ($variantseqs{$variant}[$pos-1] eq $start) {
					$variantseqs{$variant}[$pos-1] = $new;
					$variantseqs{$variant} = join("", @{$variantseqs{$variant}});
					print OUT ">$protein_id", "_", $variant, "\n", $variantseqs{$variant}, "\n";
				}
				else {
					print "failed to find residue\t$start\t$new\t$variantseqs{$variant}[$pos-1]\n";
				}
			}
		}
	}
}

