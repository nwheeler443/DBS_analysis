#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my @files = `ls Rhizosphere/*/*.faa`;		# change to Pathogenic and Environmental to complete analysis
my @done = `ls Rhizosphere/*/*.scan`;
my @doneids;
foreach my $donefile (@done) {
    if ($donefile =~ /(Rhizosphere\/.+\/.+).scan/) {
        push @doneids, $1;
    }
}

foreach my $file (@files) {
    if ($file =~ /(Rhizosphere\/.+\/.+).faa/) {
        if ($1 ~~ @doneids) {
            next;
        }
        else {
            system "hmmscan --domtblout $1.scan ../scripts/bigfiles/deltaBS.hmmlib $1.faa";
        }
    }
}