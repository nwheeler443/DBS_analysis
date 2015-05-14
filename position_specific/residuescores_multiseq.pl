#!/usr/bin/perl

## INPUT : hmm	sequence alignment

### want to be able to take an hmmalign file and an hmm file and prodice a position-wise scoring distribution for the protein pair

###### SCORING METHOD ###########

#1	2	3	4	5	6
#A	D	-	FI	C	W

#M	M	D	M	M	M
#E	E		E	E	E
#			I
#			E

#################################

# sort out last sequence problem

use warnings;
use strict;
use Data::Dumper;

my $hmmfile = shift @ARGV;
my $alignmentfile = shift @ARGV;

my %embitscores;
my %insprobabilities;
my %transprobabilities;
my $position = 0;
my $key;
my $type = "residues";

open HMM, $hmmfile;
while (<HMM>) {
	if ($_ =~ /COMPO\s+(.+)$/) {
		@{$embitscores{COMPO}} = split("  ", $1);		# choosing to use the COMPO background frequencies since im not doing bias filtering
	}
	if ($_ =~ /^\s+(\d+)\s+(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+)/) {
		$key = "pos$1";
		$position = $1;
		@{$embitscores{$key}} = split("  ", $2);
	}
	if ($_ =~ /^\s+(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+)$/) {
		$key = "pos$position";
		@{$insprobabilities{$key}} = split("  ", $1);
	}
	if ($_ =~ /^\s+(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+)$/) {
		$key = "pos$position";
		@{$transprobabilities{$key}} = split("  ", $1)
	}
}
close HMM;

#print Dumper (\@{$insprobabilities{pos95}});
#print Dumper (\@{$transprobabilities{pos359}});

my %residues = (A => 0,C => 1,D => 2,E => 3,F => 4,G => 5,H => 6,I => 7,K => 8,L => 9,M => 10,N => 11,P => 12,Q => 13,R => 14,S => 15,T => 16,V => 17,W => 18,Y => 19);

# read in alignment file, make an array called @states which has keys for insertions and matches, and put complete sequence into an array

open ALI, $alignmentfile;
my $seqnum = 1;
my %sequences;
my @states;
while (<ALI>) {
	chomp;
	if ($_ =~ /#=GC RF\s+(.+)$/) {
		push @states, split("", $1);
	}
	next if ($_ =~ /#/);
	chomp;
	if ($_ =~ /(\S+)\s+(\S+)$/) {
		push @{$sequences{$1}}, split("", $2);
	}
}
close ALI;

# go through the sequence and if a match, put it into a match sequence
# any with insertion, push that into an insertion hash with last match position as key
#record transition to that column (m->m or i->m or d->m)
# transitions = ("M", "whatever the last positions last transition was"....)

my %inserts;
my %matches;
my %transitions;

my $lastmatch = "pos0";
my $lastpoint = 0;

foreach my $key (keys(%sequences)) {
	@{$matches{$key}}[0] = "B";
}

foreach my $point (0..$#states) {
	if ($states[$point] eq "x") {
		$lastpoint = $lastpoint+1;
		$lastmatch = "pos$lastpoint";
		foreach my $key (keys(%sequences)) {
			push @{$matches{$key}}, $sequences{$key}[$point];
			if ($sequences{$key}[$point] eq "-") {
				push @{$transitions{$key}}, "D";
			}
			else {
				push @{$transitions{$key}}, "M";
			}
		}
	}
	elsif ($states[$point] eq "\.") {
		foreach my $key (keys(%sequences)) {
			if ($sequences{$key}[$point] ne ".") {
				push @{$inserts{$key}{$lastmatch}}, uc $sequences{$key}[$point];
				$transitions{$key}[$lastpoint] = "I";
			}
		}
	}
}

my %scores;

my %transitioncodes = (M => 0, I => 3, D => 5);				# m->m     m->i     m->d     i->m     i->i     d->m     d->d

# scores : transition score, match score, deletion score, insertion transition and emission score
# transition score is taken from the position before (going to a match state)

## SCORING
foreach my $key (keys(%sequences)) {
	foreach my $pos (1..$#{$matches{$key}}) {		# first position is "B"
		my $matchid = "pos$pos";
		my $prevpos = $pos-1;
		my $previd = "pos$prevpos";
		# score emission and transition
		my $residue = $matches{$key}[$pos];
		if (defined($residues{$residue})) {
			my $ref = $residues{$residue};
			push @{$scores{$key}{$matchid}}, $embitscores{$matchid}[$ref];
			push @{$scores{$key}{$matchid}}, $transprobabilities{$previd}[$transitioncodes{$transitions{$key}[$prevpos]}];
		}
		elsif ($residue eq "-") {
			push @{$scores{$key}{$matchid}}, $transprobabilities{$previd}[2];
		}
		else {
			print "residue not recognised\n";
		}
		# score inserts
		foreach my $pos (0..$#{$inserts{$key}{$matchid}}) {
			if ($pos == 0) {
				my $res = $inserts{$key}{$matchid}[$pos];
				next if ($res eq "*");
				my $insertscore = $transprobabilities{$previd}[1];
				my $insemission = $insprobabilities{$matchid}[$residues{$res}];
				push @{$scores{$key}{$matchid}}, $insertscore;
				push @{$scores{$key}{$matchid}}, $insemission;
			}
			if ($pos != 0) {
				my $res = $inserts{$key}{$matchid}[$pos];
				next if ($res eq "*");
				my $insertscore = $transprobabilities{$previd}[4];
				my $insemission = $insprobabilities{$matchid}[$residues{$res}];
				push @{$scores{$key}{$matchid}}, $insertscore;
				push @{$scores{$key}{$matchid}}, $insemission;
			}
		}
	}
	# score last position
	my $seqpos = $#{$matches{$key}};
	my $prevpos = $seqpos-1;
	my $hashid = "pos$seqpos";
	my $residue = $matches{$key}[$seqpos];
	if (defined($residues{$residue})) {
		my $ref = $residues{$residue};
		push @{$scores{$key}{$hashid}}, $embitscores{$hashid}[$ref];
	}
}

my %totalscores;

foreach my $key (keys(%sequences)) {
	foreach my $pos (1..$#{$matches{$key}}) {
		my $id = "pos$pos";
		my $total = 0;
		foreach my $score (@{$scores{$key}{$id}}) {
			$total = $total + $score;
		}
		push @{$totalscores{$key}}, $total;
	}
}

open OUT, "> scores.csv" or die "couldn't create file: $!";
foreach my $key (keys(%sequences)) {
	print OUT $key, ",";
	foreach my $pos (1..$#{$matches{$key}}) {
		print OUT $matches{$key}[$pos], ",";
	}
	print OUT "\n", $key, ",";
	foreach my $pos (0..$#{$totalscores{$key}}) {
		print OUT $totalscores{$key}[$pos], ",";
	}
	print OUT "\n";
}
close OUT;


