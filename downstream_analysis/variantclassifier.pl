#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;

my @filelist = `ls seqalignments/*.afa`;

foreach my $file (@filelist){
	open IN, $file;
	my %sequences;
	while (<IN>) {
		chomp;
		next if ($_ =~ /#/);
			if ($_ =~/(\S+)\s+(\S+)/) {
				push @{$sequences{$1}}, split("", $2);
			}
		}

	my @sequences = keys(%sequences);

	my $seq1 = $sequences[0];
	my $seq2 = $sequences[1];
	my $seq1gap = 0;
	my $seq2gap = 0;
	my $seq1del = 0;
	my $seq2del = 0;
	my $mismatches = 0;

	#print Dumper (\%sequences);

	while ($sequences{$seq1}[$#{$sequences{$seq1}}] eq "-" && $sequences{$seq2}[$#{$sequences{$seq2}}] eq "-") {
		pop $sequences{$seq1};
		pop $sequences{$seq2};
	}

	foreach my $pos (0..$#{$sequences{$seq1}}) {
		print "$sequences{$seq1}[$pos]\t$sequences{$seq2}[$pos]\n";
		if ($sequences{$seq1}[$pos] ne $sequences{$seq2}[$pos]) {
			if ($sequences{$seq1}[$pos] eq "-" | $sequences{$seq1}[$pos] eq ".") {
				$seq1gap++;
			}
			if ($sequences{$seq2}[$pos] eq "-" | $sequences{$seq2}[$pos] eq ".") {
				$seq2gap++;
			}
			else {
				$mismatches++;
				$seq1del = $seq1del + $seq1gap;
				$seq1gap = 0;
				$seq2del = $seq2del + $seq2gap;
				$seq2gap = 0;
			}
			
		}
		else {
			$seq1del = $seq1del + $seq1gap;
			$seq1gap = 0;
			$seq2del = $seq2del + $seq2gap;
			$seq2gap = 0;
		}
	}

	# truncation of seq1	truncation of seq2	deletions in seq1	deletions in seq2	mismatches
	print "$seq1gap\t$seq2gap\t$seq1del\t$seq2del\t$mismatches\n";
}