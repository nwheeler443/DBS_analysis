#!/usr/bin/perl

# make sure you align the sequences first!!!

use warnings;
use strict;
use Data::Dumper;

my $filename;
my @cutoffs = [0,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1];
my $outdir;

# send in argumants from another script to indicate input file and cutoff
foreach my $arg (@ARGV) {
	print $arg, "\t";
	if ($arg =~ /^(\S+.hits)\s+(\S+)$/) {
		$filename = $1;
		$outdir=$2;
	}
	else {print "wrong argument input\n"}
}
my %sequences;
my %gaps;
my %percentids;
my $seqname;
my $original = "";

open IN, $filename;
while(<IN>) {
    chomp;
	if ($_ =~ /=GF\ ID\ (\S+)-i\d/) {
		$original = $1;
	}
	elsif ($_ =~ /^#/) {
		next;
	}
	elsif ($_ =~ /(\S+)\s+(\S+)$/){
		my @line = split("", $2);
		foreach my $res (@line) {
			if (isAA ($res)) {
				push @{$gaps{$1}}, 0;
			}
			else {
				push @{$gaps{$1}}, 1;
			}
		}
		push @{$sequences{$1}}, @line;
	}
}

close IN;

my @nseq = keys(%sequences);
print "n=$#nseq\tl=$#{$sequences{$nseq[0]}}\n";

#print Dumper (\@{$sequences{$original}});

# create a new file for each cutoff and print original sequence to it

my %overthreshold;

# go through each hit and determine the percentage identity
open PIDS, "> troubleshooting/$original.pids.txt";
my %matchcounts;
#my $count = 0;
foreach my $seq (keys(%sequences)){
#	if ($count < 250) {
		next if ($seq eq "$original");
		my $alignmentlength = $#{$sequences{$original}} + 1;
		$matchcounts{$seq}=0;
		foreach my $residue (0..$#{$sequences{$original}}) {
			if ($gaps{$original}[$residue] == 1 && $gaps{$seq}[$residue] == 1) {
				$alignmentlength--;
			}
			elsif (uc $sequences{$original}[$residue] eq uc $sequences{$seq}[$residue]) {
				$matchcounts{$seq}++;
			}
		}
		my $percentid = $matchcounts{$seq}/$alignmentlength;
		$percentids{$seq} = $percentid;
#		if ($percentid > 0.35) {
#			$count++;
#		}
		printf PIDS "%0.2f\t$seq\n", ($matchcounts{$seq}/$alignmentlength*100);
#	}
}
close PIDS;

#produce output file once all PIDs have been calculated
foreach my $cut (@cutoffs) {
	my $outfile = "$outdir/$original.$cut.fasta";
	open OUT, "> $outfile";
	print OUT ">$original\n";
	foreach my $res (0..$#{$sequences{$original}}) {
		print OUT $sequences{$original}[$res];
	}
	print OUT "\n";
	foreach my $seq (keys(%percentids)) {
		if($percentids{$seq}>$cut) {
			print OUT ">$seq\n";
			foreach my $res (0..$#{$sequences{$seq}}) {
				print OUT $sequences{$seq}[$res];
			}
			print OUT "\n";
		}
	}
}

# remove columns with all gaps to produce aligned file
foreach my $cut (@cutoffs) {
	system "esl-reformat --mingap afa $outdir/$original.$cut.fasta > $outdir/$original.afa";
}

############################
#isAA: check if single character belongs to the IUPAC amino-acid code.
sub isAA {
    
    my $aa = shift;
    return 0 if (not defined($aa) or length($aa) != 1);
    my $iupac = 'ABCDEFGHIKLMNPQRSTVWXY';
    
    if ($aa=~/[$iupac]/i){
        return 1;
    }
	else {
		return 0;
	}
}

