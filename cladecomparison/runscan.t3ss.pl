#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my @files = `ls */*/*.faa`;
my @done = `ls */*/*.scan.t3ees`;
my @doneids;
foreach my $donefile (@done) {
    if ($donefile =~ /(.+\/.+\/.+).scan/) {
        push @doneids, $1;
    }
}

foreach my $file (@files) {
    if ($file =~ /(.+\/.+\/.+).faa/) {
        if ($1 ~~ @doneids) {
            next;
        }
        else {
            system "hmmscan --domtblout $1.scan.t3ees T3EE_Models/allT3EE.hmm $1.faa";
        }
    }
}

