#!/usr/bin/perl

### want to be able to take an hmmalign file and an hmm file and prodice a position-wise scoring distribution for the protein pair

###### SCORING METHOD ###########

#1	2	3	4	5	6
#A	D	-	FI	C	w

#M	M	D	M	M	M
#E	E		E	E	E
#			I
#			E

#################################

# sort out last sequence problem

## ONLY SCORING THE INCIDENCE OF AN INSERTION, NOT THE EMISSION ASSOCIATED WITH IT, AND NOT SCORING M->M TRANSITIONS

#my %backgroundfrequencies = (A=>0.0787945,B=>0.0151600,C=>0.0535222,D=>0.0668298,E=>0.0397062,F=>0.0695071,G=>0.0229198,0.0590092,0.0594422,0.0963728,0.0237718,0.0414386,0.0482904,0.0395639,0.0540978,0.0673417,0.0114135,0.0304133);
#{
#	f[0] = 0.0787945;		/* A */
#	f[1] = 0.0151600;		/* C */
#	f[2] = 0.0535222;		/* D */
#	f[3] = 0.0668298;		/* E */
#	f[4] = 0.0397062;		/* F */
#	f[5] = 0.0695071;		/* G */
#	f[6] = 0.0229198;		/* H */
#	f[7] = 0.0590092;		/* I */
#	f[8] = 0.0594422;		/* K */
#	f[9] = 0.0963728;		/* L */
#	f[10]= 0.0237718;		/* M */
#	f[11]= 0.0414386;		/* N */
#	f[12]= 0.0482904;		/* P */
#	f[13]= 0.0395639;		/* Q */
#	f[14]= 0.0540978;		/* R */
#	f[15]= 0.0683364;		/* S */
#	f[16]= 0.0540687;		/* T */
#	f[17]= 0.0673417;		/* V */
#	f[18]= 0.0114135;		/* W */
#	f[19]= 0.0304133;		/* Y */
#	return eslOK;
#}


#
#bg->p1    = 350./351.;
#bg->omega = 1./256.;
#bg->abc   = abc;


# logodds = eslCONST_LOG2R * log(p / bg->f[j]);			from hmmlogo



# the probability values in the hmm file are e^-(x)



use warnings;
use strict;
use Data::Dumper;

my $hmmfile = "profilehmm";
my $alignmentfile = "sequencealignment";

my %embitscores;
my %insprobabilities;
my %transprobabilities;
my $position = 0;
my $key;
my $type = "residues";

open HMM, $hmmfile;
while (<HMM>) {
	if ($_ =~ /COMPO\s+(.+)\n/) {
		@{$embitscores{COMPO}} = split("  ", $1);		# choosing to use the COMPO background frequencies since im not doing bias filtering
	}
	if ($_ =~ /^\s+(\d+)\s+(.+\d)\s+\d+\s+.\s+.$/) {
		$key = "pos$1";
		$position = $1;
		@{$embitscores{$key}} = split("  ", $2);
	}
	if ($_ =~ /^\s+(\d+.+\d+)$/) {
		$key = "pos$position";
		@{$insprobabilities{$key}} = split("  ", $1);
	}
	if ($_ =~ /^\s+(\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+)$/) {
		$key = "pos$position";
		@{$transprobabilities{$key}} = split("  ", $1)
	}
}

#print Dumper (\@{$insprobabilities{pos95}});

my %residues = (A => 0,C => 1,D => 2,E => 3,F => 4,G => 5,H => 6,I => 7,K => 8,L => 9,M => 10,N => 11,P => 12,Q => 13,R => 14,S => 15,T => 16,V => 17,W => 18,Y => 19);

# make an array called @insertions which has keys for insertions and matches

open ALI, $alignmentfile;
my @firstseq;
my @secondseq;
my @insertions;
while (<ALI>) {
	if ($_ =~ /#=GC RF\s+(.+)\n/) {
		push @insertions, split("", $1);
	}
	next if ($_ =~ /#/);
	chomp;
	if ($_ =~ /species1\s+(.+)/) {
		push @firstseq, split("", $1);
	}
	if ($_ =~ /species2\s+(.+)/) {
		push @secondseq, split("", $1);
	}
}

#print Dumper (\@insertions);

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
foreach my $point (0..$#insertions) {
	if ($insertions[$point] eq "x") {
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
	if ($insertions[$point] eq "\.") {
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

my %scores1;
my %scores2;

my %transitioncodes = (I => 3, M => 0, D => 5);				# m->m     m->i     m->d     i->m     i->i     d->m     d->d

# check that the transitions match up to the right column:
#foreach my $spot (0..$#firstmatches) {
#	print "$firsttransitions[$spot]\t$secondtransitions[$spot]\t$transitioncodes{$firsttransitions[$spot-1]}\t$transitioncodes{$secondtransitions[$spot-1]}\n"
#}

# scores : transition score, match score, deletion score, insertion transition and emission score
# transition score is taken from the position before (going to a match state)

#print "position 1, sequence 1: $firstmatches[1]\t$embitscores{pos1}[$residues{$firstmatches[1]}]\t$transprobabilities{pos0}[$transitioncodes{$firsttransitions[0]}]\n";
#print "position 1, sequence 2: $secondmatches[1]\t$embitscores{pos1}[$residues{$secondmatches[1]}]\t$transprobabilities{pos0}[$transitioncodes{$secondtransitions[0]}]\n";

## SCORING THE FIRST SEQUENCE
#score inserts at the start
#foreach my $res (@{$firstinserts{pos0}}) {
#	#print $res, "\n";
#	# need to sort out how the transition scores work here
#	my $residue = $residues{$res};
#	my $matchscore = $insprobabilities{pos0}[$residue]-$embitscores{COMPO}[$residue];
#	#print $matchscore, "\n";			# seems to be producing a mix of positive and negative values - seems promising
#	push @{$scores1{pos0}}, $matchscore;
#}
#score subsequent sequence
foreach my $pos (1..$#firstmatches) {		# first position is "B"
	my $matchid = "pos$pos";
	my $prevpos = $pos-1;
	my $previd = "pos$prevpos";
	#score transition
	my $matchscore = $transprobabilities{$previd}[$transitioncodes{$firsttransitions[$prevpos]}]+log(0.99995);			# need to get an actual null probability and work out the base for the log (just divide that by log of whatever base)
	push @{$scores1{$matchid}}, $matchscore;
	# score emission
	my $residue = $firstmatches[$pos];
	if (defined($residues{$residue})) {
		my $ref = $residues{$residue};
		push @{$scores1{$matchid}}, $embitscores{$matchid}[$ref]-$embitscores{COMPO}[$ref];		# could always add this to scores1{pos}[0] here
	}
	if ($residue eq "_") {
		push @{$scores1{$matchid}}, $transprobabilities{$matchid}[2]+log(0.99995);
	}
	foreach my $pos (0..$#{$firstinserts{$matchid}}) {
		if ($pos == 0) {
			my $res = $firstinserts{$matchid}[$pos];
			my $insertscore = $transprobabilities{$matchid}[1]+log(0.99995);
			my $insemission = $insprobabilities{$matchid}[$residues{$res}]-$embitscores{COMPO}[$residues{$res}];
			push @{$scores1{$matchid}}, $insertscore;
			push @{$scores1{$matchid}}, $insemission;
		}
		if ($pos != 0) {
			my $res = $firstinserts{$matchid}[$pos];
			my $insertscore = $transprobabilities{$matchid}[4]+log(0.99995);
			my $insemission = $insprobabilities{$matchid}[$residues{$res}]-$embitscores{COMPO}[$residues{$res}];
			push @{$scores1{$matchid}}, $insertscore;
			push @{$scores1{$matchid}}, $insemission;
		}
	}
}
my $seqpos = $#firstmatches;
my $prevpos = $#firstmatches-1;
my $hashid = "pos$seqpos";
#my $matchscore = $transbitscores{$hashid}[$transitioncodes{$firsttransitions[$pos-1]}];
#push @{$scores1{$hashid}}, $matchscore;
my $residue = $firstmatches[$seqpos];
if (defined($residues{$residue})) {
	my $ref = $residues{$residue};
	push @{$scores1{$hashid}}, $embitscores{$hashid}[$ref];
}

#print Dumper (\%scores1);

## SCORING THE SECOND SEQUENCE
#score inserts at the start
#foreach my $res (@{$secondinserts{pos0}}) {
#	#print $res, "\n";
#	# need to sort out how the transition scores work here
#	my $residue = $residues{$res};
#	my $matchscore = $insprobabilities{pos0}[$residue]-$embitscores{COMPO}[$residue];
#	#print $matchscore, "\n";			# seems to be producing a mix of positive and negative values - seems promising
#	push @{$scores2{pos0}}, $matchscore;
#}
foreach my $pos (1..$#secondmatches) {
	my $matchid = "pos$pos";
	my $prevpos = $pos-1;
	my $previd = "pos$prevpos";
#	#score transition
	my $matchscore = $transprobabilities{$previd}[$transitioncodes{$secondtransitions[$pos]}]+log(0.99995);
	push @{$scores2{$matchid}}, $matchscore;
#	# score emission
	my $residue = $secondmatches[$pos];
	if (defined($residues{$residue})) {
		my $ref = $residues{$residue};
		push @{$scores2{$matchid}}, $embitscores{$matchid}[$ref]-$embitscores{COMPO}[$ref];
	}
	if ($residue eq "_") {
		push @{$scores2{$matchid}}, $transprobabilities{$matchid}[2]+log(0.99995);
	}
	# score any inserts
	foreach my $pos (0..$#{$secondinserts{$matchid}}) {
		my $res = $secondinserts{$matchid}[$pos];
		my $insertscore = $transprobabilities{$matchid}[1]+log(0.99995);
		my $insemission = $insprobabilities{$matchid}[$residues{$res}]-$embitscores{COMPO}[$residues{$res}];
		#print $res, "\t", $hashid, "\t", $insertscore, "\t", $insemission, "\n";
		push @{$scores2{$matchid}}, $insertscore;
		push @{$scores2{$matchid}}, $insemission;
	}
}
$seqpos = $#secondmatches;
$prevpos = $#secondmatches-1;
$hashid = "pos$seqpos";
#my $matchscore = $transbitscores{$hashid}[$transitioncodes{$secondtransitions[$pos-1]}];
#push @{$scores2{$hashid}}, $matchscore;
$residue = $firstmatches[$seqpos];
if (defined($residues{$residue})) {
	my $ref = $residues{$residue};
	push @{$scores2{$hashid}}, $embitscores{$hashid}[$ref];
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
	push @totalscores1, ($total*(-1));
	#print "\n";
}
foreach my $pos (1..$#secondmatches) {
	$key = "pos$pos";
	my $total = 0;
	foreach my $score (@{$scores2{$key}}) {
		$total = $total + $score;
	}
	push @totalscores2, ($total*(-1));
}

#print Dumper (\@totalscores1);
#print Dumper (\@totalscores2);

open OUT, "> scores.csv" or die "couldn't create file: $!";
foreach my $pos (1..$#firstmatches) {
	print OUT $firstmatches[$pos], ",";
}
print OUT "\n";
foreach my $pos (0..$#totalscores1) {
	print OUT $totalscores1[$pos], ",";
}
print OUT "\n";
foreach my $pos (1..$#secondmatches) {
	print OUT $secondmatches[$pos], ",";
}
print OUT "\n";
foreach my $pos (0..$#totalscores2) {
	print OUT $totalscores2[$pos], ",";
}

close OUT;

#print Dumper (\@{$scores1{pos2}});
#print Dumper (\@{$scores2{pos2}});

#print "$#totalscores1\t$#totalscores2";


