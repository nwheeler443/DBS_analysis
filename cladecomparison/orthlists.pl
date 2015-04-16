#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my @pathogenic = `ls Pathogenic/*/*.scan`;
my @environmental = `ls Environmental/*/*.scan`;
my @rhizosphere = `ls Rhizosphere/*/*.scan`;

foreach my $table (@pathogenic) {
	if ($table =~ /(Pathogenic\/.+\/.+).scan/) {
		if ($1 eq "Pathogenic/Pto_DC3000_P/Pto_DC3000_P_AE16853") {
			next;
		}
		else {
			if (-e "$1.orths") {
				next;
			}
			else {
				system "~/Dropbox/scripts/orthlist.pl Pathogenic/Pto_DC3000_P/Pto_DC3000_P_AE16853.faa $1.faa Pathogenic/Pto_DC3000_P/Pto_DC3000_P_AE16853.scan $1.scan $1.orths";
			}
		}
	}
}

foreach my $table (@environmental) {
	if ($table =~ /(Environmental\/.+\/.+).scan/) {
		if ($1 eq "Environmental/Pfl_Pf0-1_E/Pfl_Pf0-1_E_NC_007492") {
			next;
		}
		else {
			if (-e "$1.orths") {
				next;
			}
			else {
				system "~/Dropbox/scripts/orthlist.pl Environmental/Pfl_Pf0-1_E/Pfl_Pf0-1_E_NC_007492.faa $1.faa Environmental/Pfl_Pf0-1_E/Pfl_Pf0-1_E_NC_007492.scan $1.scan $1.orths";
			}
		}
	}
}

foreach my $table (@rhizosphere) {
	if ($table =~ /(Rhizosphere\/.+\/.+).scan/) {
		if ($1 eq "Rhizosphere/Pfl_PCL1751_R/Pfl_PCL1751_R_CP010896") {
			next;
		}
		else {
			if (-e "$1.orths") {
				next;
			}
			else {
				system "~/Dropbox/scripts/orthlist.pl Rhizosphere/Pfl_PCL1751_R/Pfl_PCL1751_R_CP010896.faa $1.faa Rhizosphere/Pfl_PCL1751_R/Pfl_PCL1751_R_CP010896.scan $1.scan $1.orths";
			}
		}
	}
}