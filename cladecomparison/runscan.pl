#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my @files = `ls Environmental/*/*.faa`;
my @done = `ls Environmental/*/*.scan`;
my @doneids;
foreach my $donefile (@done) {
    if ($donefile =~ /(Environmental\/.+\/.+).scan/) {
        push @doneids, $1;
    }
}

foreach my $file (@files) {
    if ($file =~ /(Environmental\/.+\/.+).faa/) {
        if ($1 ~~ @doneids) {
            next;
        }
        else {
            system "hmmscan --domtblout $1.scan ../scripts/bigfiles/deltaBS.hmmlib $1.faa";
        }
    }
}