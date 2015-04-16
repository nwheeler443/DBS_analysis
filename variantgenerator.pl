#!/usr/bin/perl

# given the humsavar database and a FASTA file of all the original sequences in the database, will genrate variant sequences for all of them

use warnings;
use strict;
use Data::Dumper;

my %variants;
open IN, "uniprot.humsavar.2011_09.provean.scored";
while (<IN>) {
	if ($_ =~ /^(\S+)\s+(\S+)/) {
		push @{$variants{$1}}, $2;
	}
}

open OUT, "> query_sequences.txt";

my @sequences = `ls fastafiles/A*.fasta`;
foreach my $file (@sequences) {
	if ($file =~ /fastafiles\/(.+).fasta/) {
		my $protein_id = $1;
		my @sequence;
		open FASTA, "fastafiles/$1.fasta";			# fasta file
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

@sequences = `ls fastafiles/P*.fasta`;
foreach my $file (@sequences) {
	if ($file =~ /fastafiles\/(.+).fasta/) {
		my $protein_id = $1;
		my @sequence;
		open FASTA, "fastafiles/$1.fasta";			# fasta file
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
				print "failed to find residue\t$start\t$new\t$variantseqs{$variant}[$pos-1]\t$protein_id\n";
			}
		}
	}
}

foreach my $number (0..9) {
	my @sequences = `ls fastafiles/Q$number*.fasta`;
	foreach my $file (@sequences) {
		if ($file =~ /fastafiles\/(.+).fasta/) {
			my $protein_id = $1;
			my @sequence;
			open FASTA, "fastafiles/$1.fasta";			# fasta file
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
