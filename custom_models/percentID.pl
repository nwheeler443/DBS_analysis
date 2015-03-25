#!/usr/bin/perl

# make sure you align the sequences first!!!
# have a sequence at the top labelled original

use warnings;
use strict;
use Data::Dumper;

my %sequences;
my $seqname;
open IN, shift @ARGV;
while(<IN>) {
    chomp;
    if($_ =~/>(\S+)/){
        $seqname=$1;
    }
    else {push @{$sequences{$seqname}}, split("", $_)}
}
close IN;

#print Dumper (\%sequences);

my @cutoffs= @ARGV;

foreach my $cut (@cutoffs) {
    my $outfile = "$cut.fasta";
    open OUT, ">", $outfile;
    print OUT ">original\n";
    foreach my $res (0..$#{$sequences{original}}) {
        if (isAA($sequences{original}[$res])) {
            print OUT $sequences{original}[$res];
        }
    }
    print OUT "\n";
    close OUT;
}

my %overthreshold;

open PIDS, "> pids.txt";
my %matchcounts;
foreach my $seq (keys(%sequences)){
    next if ($seq eq "original");
    my $alignmentlength = $#{$sequences{original}} + 1;
    $matchcounts{$seq}=0;
    foreach my $residue (0..$#{$sequences{original}}) {
#        if ($sequences{original}[$residue] eq "-" && $sequences{$seq}[$residue] eq "-") {
#            $alignmentlength--;
#        }
#        elsif (uc $sequences{original}[$residue] eq uc $sequences{$seq}[$residue]) {
#            $matchcounts{$seq}++;
#        }
        if (isAA ($sequences{original}[$residue]) && isAA ($sequences{$seq}[$residue]) && uc $sequences{original}[$residue] eq uc $sequences{$seq}[$residue]) {
            $matchcounts{$seq}++;
        }
        elsif ( !isAA ($sequences{original}[$residue]) && !isAA($sequences{$seq}[$residue]) ) {
            $alignmentlength--;
        }
        
    }
    printf PIDS "%0.2f\t$seq\n", ($matchcounts{$seq}/$alignmentlength*100);
    foreach my $cut (@cutoffs) {
        my $outfile = "$cut.fasta";
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

