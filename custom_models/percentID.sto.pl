#!/usr/bin/perl

# TO DO:
# get rid of the original naming requirement
# re-format read-in of alignment

# make sure you align the sequences first!!!
# have a sequence at the top labelled original

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

foreach my $arg (@ARGV) {
	print $arg, "\n";
	if ($arg =~ /^(jackhmmerresults\/\S+.hits)\n.(\S+)$/) {
		$filename = $1;
		@cutoffs = $2;
	}
	else {print "wrong\n"}
}

my %sequences;
my $seqname;
my $original = "";

open IN, $filename;
while(<IN>) {
    chomp;
	if ($_ =~ /=GF\ ID\ (\S+)-i\d/) {
		$original = $1;
	}
	elsif ($_ =~ /#/) {
		next;
	}
else {
	if ($_ =~ /(\S+)\s+(\S+)$/){
		push @{$sequences{$1}}, split("", $2);
	}
}
}
close IN;

#print Dumper (\%sequences);

foreach my $cut (@cutoffs) {
    my $outfile = "filtered/$original.$cut.fasta";
    open OUT, ">", $outfile;
    print OUT ">$original\n";
    foreach my $res (0..$#{$sequences{$original}}) {
        if (isAA($sequences{$original}[$res])) {
            print OUT $sequences{$original}[$res];
        }
    }
    print OUT "\n";
    close OUT;
}

my %overthreshold;

open PIDS, "> troubleshooting/$original.pids.txt";
my %matchcounts;
foreach my $seq (keys(%sequences)){
    next if ($seq eq "$original");
    my $alignmentlength = $#{$sequences{$original}} + 1;
    $matchcounts{$seq}=0;
    foreach my $residue (0..$#{$sequences{$original}}) {
        if (isAA ($sequences{$original}[$residue]) && isAA ($sequences{$seq}[$residue]) && uc $sequences{$original}[$residue] eq uc $sequences{$seq}[$residue]) {
            $matchcounts{$seq}++;
        }
        elsif ( !isAA ($sequences{$original}[$residue]) && !isAA($sequences{$seq}[$residue]) ) {
            $alignmentlength--;
        }
        
    }
	#print $alignmentlength;
    printf PIDS "%0.2f\t$seq\n", ($matchcounts{$seq}/$alignmentlength*100);
    foreach my $cut (@cutoffs) {
        my $outfile = "filtered/$original.$cut.fasta";
        open OUT, ">>", $outfile;
        if(($matchcounts{$seq}/$alignmentlength)>$cut) {
            print OUT ">$seq\n";
            foreach my $res (0..$#{$sequences{$seq}}) {
                if (isAA ($sequences{$seq}[$res])) {
                    print OUT $sequences{$seq}[$res];
                }
            }
           print OUT "\n";
        }
        close OUT;
    }
}
close PIDS;

############################
#isAA: check if single character belongs to the IUPAC amino-acid code.
sub isAA {
    
    my $aa = shift;
    return 0 if (not defined($aa) or length($aa) != 1);
    my $iupac = 'ABCDEFGHIKLMNPQRSTVWXYZ';
    
    if ($aa=~/[$iupac]/i){
        return 1;
    }
    
    return 0;
}

