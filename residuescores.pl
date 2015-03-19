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
my $name1;
my $name2;
my @firstseq;
my @secondseq;
my @states;
while (<ALI>) {
	chomp;
	if ($_ =~ /#=GC RF\s+(.+)$/) {
		push @states, split("", $1);
	}
	next if ($_ =~ /#/);
	chomp;
	if ($_ =~ /(\S+)\s+(\S+)$/) {
		if ($seqnum == 1) {
			push @firstseq, split("", $2);
			$seqnum = 2;
			$name1 = $1;
		}
		elsif ($seqnum == 2) {
			push @secondseq, split("", $2);
			$name2 = $1;
			$seqnum=1;
		}
	}
}
close ALI;

#print Dumper (@states);
#print Dumper (@firstseq);
#print Dumper (\@secondseq);

# go through the sequence and if a match, put it into a match sequence
# any with insertion, push that into an insertion hash with last match position as key
#record transition to that column (m->m or i->m or d->m)
# transitions = ("M", "whatever the last positions last transition was"....)

my %firstinserts;
my %secondinserts;
my $lastmatch = "pos0";
my $lastpoint = 0;
my @firstmatches = ("B");
my @secondmatches = ("B");
my @firsttransitions;
my @secondtransitions;
foreach my $point (0..$#states) {
	if ($states[$point] eq "x") {
		$lastpoint = $lastpoint+1;
		my $seqpos = $point + 1;
		$lastmatch = "pos$seqpos";
		push @firstmatches, $firstseq[$point];
		push @secondmatches, $secondseq[$point];
		if ($firstseq[$point] eq "-") {
			push @firsttransitions, "D";
		}
		else {
			push @firsttransitions, "M";
		}
		if ($secondseq[$point] eq "-") {
			push @secondtransitions, "D";
		}
		else {
			push @secondtransitions, "M";
		}
	}
	if ($states[$point] eq "\.") {
		if ($firstseq[$point] ne "\.") {
			push @{$firstinserts{$lastmatch}}, uc $firstseq[$point];
			$firsttransitions[$lastpoint] = "I";			#want to over-write the state for the last point as ending in an insertion
		}
		if ($secondseq[$point] ne "\.") {
			push @{$secondinserts{$lastmatch}}, uc $secondseq[$point];
			$secondtransitions[$lastpoint] = "I";
		}
	}
}

#print Dumper (\@firsttransitions);
#print Dumper (\@secondtransitions);
#print Dumper (\%firstinserts);
#print Dumper (\%secondinserts);
#print Dumper (\@firstmatches);
#
my %scores1;
my %scores2;

my %transitioncodes = (M => 0, I => 3, D => 5);				# m->m     m->i     m->d     i->m     i->i     d->m     d->d

# check that the transitions match up to the right column:
#foreach my $spot (0..$#firstmatches) {
#	print "$firsttransitions[$spot]\t$secondtransitions[$spot]\t$transitioncodes{$firsttransitions[$spot-1]}\t$transitioncodes{$secondtransitions[$spot-1]}\n"
#}

# scores : transition score, match score, deletion score, insertion transition and emission score
# transition score is taken from the position before (going to a match state)

## SCORING THE FIRST SEQUENCE
foreach my $pos (1..$#firstmatches) {		# first position is "B"
	my $matchid = "pos$pos";
	my $prevpos = $pos-1;
	my $previd = "pos$prevpos";
	# score emission and transition
	my $residue = $firstmatches[$pos];
	if (defined($residues{$residue})) {
		my $ref = $residues{$residue};
		#print "$matchid\t", $embitscores{COMPO}[$ref]-$embitscores{$matchid}[$ref], "\t", log(0.99995)/log(2)-$transprobabilities{$previd}[$transitioncodes{$firsttransitions[$prevpos]}], "\n";
		push @{$scores1{$matchid}}, $embitscores{COMPO}[$ref]-$embitscores{$matchid}[$ref];
		push @{$scores1{$matchid}}, log(0.99995)/log(2)-$transprobabilities{$previd}[$transitioncodes{$firsttransitions[$prevpos]}];			# need to get an actual null probability and work out the base for the log (just divide that by log of whatever base)
	}
	elsif ($residue eq "-") {
		push @{$scores1{$matchid}}, log(0.99995)/log(2)-$transprobabilities{$previd}[2];
	}
	else {
		print "residue not recognised\n";
	}
	# score inserts
	my $runningtotal = 0;
	foreach my $pos (0..$#{$firstinserts{$matchid}}) {
		if ($pos == 0) {
			my $res = $firstinserts{$matchid}[$pos];
			my $insertscore = log(0.99995)/log(2)-$transprobabilities{$previd}[1];
			my $insemission = $embitscores{COMPO}[$residues{$res}]-$insprobabilities{$matchid}[$residues{$res}];
			push @{$scores1{$matchid}}, $insertscore;
			push @{$scores1{$matchid}}, $insemission;
			$runningtotal = $runningtotal+$insertscore+$insemission;
		}
		if ($pos != 0) {
			my $res = $firstinserts{$matchid}[$pos];
			my $insertscore = log(0.99995)/log(2)-$transprobabilities{$previd}[4];
			my $insemission = $embitscores{COMPO}[$residues{$res}]-$insprobabilities{$matchid}[$residues{$res}];
			push @{$scores1{$matchid}}, $insertscore;
			push @{$scores1{$matchid}}, $insemission;
			$runningtotal = $runningtotal+$insertscore+$insemission;
		}
	}
}
# score last position
my $seqpos = $#firstmatches;
my $prevpos = $#firstmatches-1;
my $hashid = "pos$seqpos";
my $residue = $firstmatches[$seqpos];
if (defined($residues{$residue})) {
	my $ref = $residues{$residue};
	push @{$scores1{$hashid}}, $embitscores{COMPO}[$residues{$ref}]-$embitscores{$hashid}[$ref];
}

#print Dumper (\%scores1);

# SCORING THE SECOND SEQUENCE
foreach my $pos (1..$#secondmatches) {
	my $matchid = "pos$pos";
	my $prevpos = $pos-1;
	my $previd = "pos$prevpos";
	# score emission and transition
	my $residue = $secondmatches[$pos];
	if (defined($residues{$residue})) {
		my $ref = $residues{$residue};
		push @{$scores2{$matchid}}, $embitscores{COMPO}[$ref]-$embitscores{$matchid}[$ref];
		push @{$scores2{$matchid}}, log(0.99995)/log(2)-$transprobabilities{$previd}[$transitioncodes{$secondtransitions[$prevpos]}];
	}
	elsif ($residue eq "-") {
		push @{$scores2{$matchid}}, log(0.99995)/log(2)-$transprobabilities{$previd}[2];
	}
	else {
		print "residue not recognised\n";
	}
	# score any inserts
	foreach my $pos (0..$#{$secondinserts{$matchid}}) {
		if ($pos == 0) {
			my $res = $secondinserts{$matchid}[$pos];
			my $insertscore = log(0.99995)/log(2)-$transprobabilities{$previd}[1];
			my $insemission = $embitscores{COMPO}[$residues{$res}]-$insprobabilities{$matchid}[$residues{$res}];
			push @{$scores2{$matchid}}, $insertscore;
			push @{$scores2{$matchid}}, $insemission;
		}
		if ($pos != 0) {
			my $res = $secondinserts{$matchid}[$pos];
			my $insertscore = log(0.99995)/log(2)-$transprobabilities{$previd}[4];			# insert extension
			my $insemission = $embitscores{COMPO}[$residues{$res}]-$insprobabilities{$matchid}[$residues{$res}];
			push @{$scores2{$matchid}}, $insertscore;
			push @{$scores2{$matchid}}, $insemission;
		}
	}
}
$seqpos = $#secondmatches;
$prevpos = $#secondmatches-1;
$hashid = "pos$seqpos";
$residue = $firstmatches[$seqpos];
if (defined($residues{$residue})) {
	my $ref = $residues{$residue};
	push @{$scores2{$hashid}}, $embitscores{COMPO}[$residues{$ref}]-$embitscores{$hashid}[$ref];
}

#print Dumper (\%scores2);

my @totalscores1;
my @totalscores2;

# don't know what this step is doing - hope its summing emission and transition scores

foreach my $pos (1..$#firstmatches) {
	$key = "pos$pos";
	#print $key, "\t";
	my $total = 0;
	foreach my $score (@{$scores1{$key}}) {
		#print $score, "\t";
		$total = $total + $score;
	}
	push @totalscores1, $total;
	#print "\n";
}
foreach my $pos (1..$#secondmatches) {
	$key = "pos$pos";
	my $total = 0;
	foreach my $score (@{$scores2{$key}}) {
		$total = $total + $score;
	}
	push @totalscores2, $total;
}

#print Dumper (\@totalscores1);
#print Dumper (\@totalscores2);

open OUT, "> scores.csv" or die "couldn't create file: $!";
print OUT $name1, ",";
foreach my $pos (1..$#firstmatches) {
	print OUT $firstmatches[$pos], ",";
}
print OUT "\n,";
foreach my $pos (0..$#totalscores1) {
	print OUT $totalscores1[$pos], ",";
}
print OUT "\n", $name2, ",";
foreach my $pos (1..$#secondmatches) {
	print OUT $secondmatches[$pos], ",";
}
print OUT "\n,";
foreach my $pos (0..$#totalscores2) {
	print OUT $totalscores2[$pos], ",";
}

close OUT;

#print Dumper (\@{$scores1{pos2}});
#print Dumper (\@{$scores2{pos2}});

#print "$#totalscores1\t$#totalscores2";


