#!/usr/bin/perl

## you would run this from Terminal by typing: ./fetchFASTAinfo.pl fastafile.fasta

use warnings;
use strict;

## this says take the stuff I've written into the Terminal window after the name of the program and put it an array called @ARGV, then take the first item in this array (that's what shift is for)

my $fileName = shift @ARGV;

##

open(F, "< $fileName");
while (my $f = <F>){
    chomp($f);
    if ($f=~/\>(\w+)/){
		print $1, "\t", $1, "\t", $1, "\n";
	}
    
}
close(F);

exit(0); 



