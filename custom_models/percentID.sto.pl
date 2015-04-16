#!/usr/bin/perl

# make sure you align the sequences first!!!

use warnings;
use strict;
use Data::Dumper;

#my $filename = shift @ARGV;
#my @cutoffs = @ARGV;
#
#print $filename;
#print Dumper @cutoffs;

my $filename;
my @cutoffs;

# send in argumants from another script to indicate input file and cutoff
foreach my $arg (@ARGV) {
	print $arg, "\n";
	if ($arg =~ /^(\S+\/\S+.hits)\s+(\S+)$/) {
		$filename = $1;
		@cutoffs = $2;
	}
	else {print "wrong argument input\n"}
}

my %sequences;
my %gaps;
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

#print Dumper (\%sequences);

# create a new file for each cutoff and print original sequence to it
foreach my $cut (@cutoffs) {
    my $outfile = "filtered/$original.$cut.fasta";
    open OUT, ">", $outfile;
    print OUT ">$original\n";
    foreach my $res (0..$#{$sequences{$original}}) {
        print OUT $sequences{$original}[$res];
    }
    print OUT "\n";
    close OUT;
}

my %overthreshold;

# go through each hit and determine the percentage identity
open PIDS, "> troubleshooting/$original.pids.txt";
my %matchcounts;
foreach my $seq (keys(%sequences)){
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
    printf PIDS "%0.2f\t$seq\n", ($matchcounts{$seq}/$alignmentlength*100);
    foreach my $cut (@cutoffs) {
        my $outfile = "filtered/$original.$cut.fasta";
        open OUT, ">>", $outfile;
        if($percentid>$cut) {
            print OUT ">$seq\n";
            foreach my $res (0..$#{$sequences{$seq}}) {
				print OUT $sequences{$seq}[$res];
            }
           print OUT "\n";
        }
        close OUT;
    }
}
close PIDS;

# remove columns with all gaps to produce aligned file
foreach my $cut (@cutoffs) {
	system "esl-reformat --mingap afa filtered/$original.$cut.fasta > filtered/$original.$cut.fasta";
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
    
    return 0;
}

