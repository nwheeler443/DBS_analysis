#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# input : ~Dropbox/scripts/cladecomparison/runscan.pl group1 group2 group3

my @groups = @ARGV;

my @files;
my @done;

foreach my $group (@groups) {
	push @files, `ls $group/*.faa`;
	push @done, `ls $group/*.scan`;
}

my @doneids;
foreach my $donefile (@done) {
    if ($donefile =~ /(.+\/.+).scan/) {
        push @doneids, $1;
    }
}

foreach my $file (@files) {
    if ($file =~ /(.+\/.+).faa/) {
        if ($1 ~~ @doneids) {
            next;
        }
        else {
            system "hmmscan --domtblout $1.scan ~/Dropbox/scripts/bigfiles/deltaBS.hmmlib $1.faa";
        }
    }
}