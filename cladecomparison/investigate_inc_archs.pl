#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

my $file = shift @ARGV;
my %archs;

open IN, "$file";
while (<IN>) {
	chomp;
	my @split = split(";", $_);
	@{$archs{$split[0]}} = @split[2..3];
}
close IN;

open OUT, "> inc_arch stats";
print OUT "Matching\tUnique to gene 1\tUnique to gene 2\n";

foreach my $orth (keys(%archs)) {
	my @arch1 = split(",", ${$archs{$orth}}[0]);
	my @arch2 = split(",", ${$archs{$orth}}[1]);
	my $matchcount = 0;
#	print Dumper (\@arch1);
#	print Dumper (\@arch2);
	foreach my $query (@arch1) {
		my $matched = 0;
		foreach my $dom (0..$#arch2) {
			if ($matched == 0) {
				#print "$query\t$arch2[$dom-$matchcount]\n";
				if ($query eq $arch2[$dom-$matchcount]) {
					#print "match\t$matchcount\n";
					#				print Dumper (\@arch2);
					splice (@arch2, $dom-$matchcount, 1);				# want to double check this works
					#				print Dumper (\@arch2);
					$matchcount++;
					$matched = 1;
			}
			}
		}
	}
	my $unique1 = $#arch1 - $matchcount + 1;
	my $unique2 = $#arch2 + 1;
	my $matching = $matchcount;
	print OUT "$matching\t$unique1\t$unique2\n";
}

