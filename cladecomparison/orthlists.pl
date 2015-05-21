#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;

my @group;
my @rep;

GetOptions ("group=s" => \@group,
"rep=s" => \@rep);

if (not @group) {
	print "please specify groups for comparison using --group\n";
}
if (not @rep) {
	print "please specify reference genomes using --rep\n";
}

foreach my $num (0..$#group) {
	my $group = $group[$num];
	my $rep = $rep[$num];
	my @members = `ls $group/*.faa`;
	my $comp;
	foreach my $member (@members) {
		if ($member =~ /$group\/(.+).faa/) {
			$comp = $1;
		}
		if ($member =~ /$rep/) {
			if ($num > 0) {
				system "~/Dropbox/scripts/orthlist.pl $group[0]/$rep[0].faa $group/$rep.faa $group[0]/$rep[0].scan $group/$rep.scan $group/reporths.orths";
			}
			else {
				next;
			}
		}
		else {
			system "~/Dropbox/scripts/orthlist.pl $group/$rep.faa $group/'$comp'.faa $group/$rep.scan $group/'$comp'.scan $group/'$comp'.orths";
		}
	}
}