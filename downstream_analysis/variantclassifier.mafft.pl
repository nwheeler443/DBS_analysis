#!/usr/bin/perl

use warnings;
use strict;

use Data::Dumper;

## forgot to add in alternative starts

my @filelist = `ls seqalignments/*.afa`;

open OUT, "> variantclassifications.txt";
print OUT "species1\tspecies2\talternate start in seq1\talternate start in seq2\ttruncation of seq1\ttruncation of seq2\tdeletions in seq1\tdeletions in seq2\tmismatches\n";

my $seqname;
foreach my $file (@filelist){
	open IN, $file;
	my %sequences;
	my @sequences;
	while (<IN>) {
		chomp;
		if ($_ =~/>(\S+)/) {
			$seqname=$1;
			push @sequences, $1;
		}
		else {
			push @{$sequences{$seqname}}, split("", $_);
		}
		}
	
#	print Dumper (\@sequences);

	my $seq1 = $sequences[0];
	my $seq2 = $sequences[1];
	my $altstart = "Y";
	my $seq1altstart = 0;
	my $seq2altstart = 0;
	my $seq1gap = 0;
	my $seq2gap = 0;
	my $seq1del = 0;
	my $seq2del = 0;
	my $mismatches = 0;

#	print Dumper (\%sequences);

	while ($sequences{$seq1}[$#{$sequences{$seq1}}] eq "-" && $sequences{$seq2}[$#{$sequences{$seq2}}] eq "-") {
		pop $sequences{$seq1};
		pop $sequences{$seq2};
	}

	foreach my $pos (0..$#{$sequences{$seq1}}) {
		## find alternative starts
#		print "$sequences{$seq1}[$pos]\t$sequences{$seq2}[$pos]\n";
		if ($sequences{$seq1}[$pos] ne $sequences{$seq2}[$pos]) {
			if ($sequences{$seq1}[$pos] eq "-" | $sequences{$seq1}[$pos] eq ".") {
				if ($altstart eq "Y") {
					$seq1altstart++;
				}
				else {
					$seq1gap++;
				}
			}
			if ($sequences{$seq2}[$pos] eq "-" | $sequences{$seq2}[$pos] eq ".") {
				if ($altstart eq "Y") {
					$seq2altstart++;
				}
				else {
					$seq2gap++;
				}
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
			$altstart = "N";
		}
	}

	# truncation of seq1	truncation of seq2	deletions in seq1	deletions in seq2	mismatches
	print OUT "$seq1\t$seq2\t$seq1altstart\t$seq2altstart\t$seq1gap\t$seq2gap\t$seq1del\t$seq2del\t$mismatches\n";
}