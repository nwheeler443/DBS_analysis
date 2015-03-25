#!/usr/bin/perl 

## INPUT: fasta1 fasta2 scan1 scan2 outfile

use warnings;
use strict; 

my $fasta1 = shift @ARGV;
my $fasta2 = shift @ARGV;

my $pfamannot1 = shift @ARGV;
#my $pfamannot1 = "$fasta1-pfam_hmmscan1.tbl";
my $pfamannot2 = shift @ARGV;
#my $pfamannot2 = "$fasta2-pfam_hmmscan1.tbl";

system "head $fasta1";
system "head $fasta2";
system "head $pfamannot1";
system "head $pfamannot2";

my $hmm_lib_path = "~/Dropbox/scripts";
my $tmp_dir = "/tmp";

if(not defined($pfamannot1)){
	$pfamannot1 = "$fasta1-pfam_hmmscan1.tbl";
	print STDERR "Running hmmscan on [$fasta1] sequences with Pfam HMMs...\n";
	system("hmmscan -o /dev/null --noali --domtblout $pfamannot1 --cut_tc $hmm_lib_path/deltaBS.hmmlib $fasta1 1>&2") == 0 or die "hmmscan failed: $!";
}
else{
	print STDERR"skipping hmmscan\'ing [$fasta1]. Using [$pfamannot1] instead.\n";
}

if(not defined($pfamannot2)){
	$pfamannot2 = "$fasta2-pfam_hmmscan1.tbl";
	print STDERR "Running hmmscan on [$fasta2] sequences with Pfam HMMs...\n";
	system("hmmscan -o /dev/null --noali --domtblout $pfamannot2 --cut_tc $hmm_lib_path/deltaBS.hmmlib $fasta2 1>&2") == 0 or die "hmmscan failed: $!";
}
else{
	print STDERR "skipping hmmscan\'ing [$fasta2]. Using [$pfamannot2] instead.\n";
}

my $scan1 = &parse_hmmscan_tbl($pfamannot1);
my $scan2 = &parse_hmmscan_tbl($pfamannot2);

sub parse_hmmscan_tbl {
	my ($tbl) = @_;
	my %hits;
	open TBL, "<", $tbl;
	while(<TBL>){
		next if($_ =~ /^#/);
			chomp;
		my @splat = split /\s+/;
		#domain, eval, score, start, end
		push @{$hits{$splat[3]}}, [$splat[1], $splat[6], $splat[7], $splat[17], $splat[18]];
	}
	close TBL;
	return \%hits;
}

my %orths;
print STDERR "Since no ortholog list provided, predicting orthologs with phmmer...\n";
my $ref;
$ref = &predict_orths_phmmer($fasta1, $fasta2);
%orths = %$ref;
open OUT, ">", shift @ARGV;
foreach my $key (keys(%orths)){
	print OUT $key,"\t",$orths{$key},"\n";
}
print STDERR "done. Ad hoc ortholog list printed to orthlist.dbs\n";
close OUT;

sub predict_orths_phmmer {
	my($fasta1, $fasta2) = ($_[0], $_[1]);
	my ($line, $seq, $name);
	my ($f1, $f2, $l1, $l2);
	my $f1_tmp = "$tmp_dir/DELTABS_f1.fa";
	my $f2_tmp = "$tmp_dir/DELTABS_f2.fa";
	my $f1_tbl = "$tmp_dir/DELTABS_fa1.tbl";
	my $f2_tbl = "$tmp_dir/DELTABS_fa2.tbl";
	my $f1_domtbl = "$tmp_dir/DELTABS_fa1.domtbl";
	my $f2_domtbl = "$tmp_dir/DELTABS_fa2.domtbl";
	my %orths;
	my %pairs;
	my %coverage;
	
	#read in fasta files
	#$f1 = &read_fasta($fasta1);
	#$f2 = &read_fasta($fasta2);
	$l1 = &find_seq_lengths($fasta1);
	$l2 = &find_seq_lengths($fasta2);
	
	my $l1_no = scalar keys(%$l1);
	
	if (defined($_[2]) && defined($_[3])){
		($f1_domtbl, $f2_domtbl)=($_[2], $_[3]);
	}
	else {
		#loop through f1 searching each seq against fasta2
		#reciprocally search seqs from f2 which are hit at
		#E < 0.1
		system("phmmer -o $f1_tbl.phmmer --noali --domtblout $f1_domtbl --tblout $f1_tbl -E 0.000001 $fasta1 $fasta2") == 0 or die "phmmer failed: $!";
		system("phmmer -o $f2_tbl.phmmer --noali --domtblout $f2_domtbl --tblout $f2_tbl -E 0.000001 $fasta2 $fasta1") == 0 or die "phmmer failed: $!";
	}
	
	#	my $f1_hits = &parse_phmmer_tbl2($f1_tbl);
	#	my $f2_hits = &parse_phmmer_tbl2($f2_tbl);
	my ($f1_hits, $f1_lengths) = &parse_phmmer_tbl3($f1_domtbl);
	my ($f2_hits, $f2_lengths) = &parse_phmmer_tbl3($f2_domtbl);
	
	foreach my $id1 (keys(%$l1)){
		
		my $max = 1;
		#print "predict_orths_phmmer:[$id1]\t";
		#print "\t2:1[$f2_hits->{$id1}]\n";
		
		next if(!$f2_hits->{$id1});
		my $f2_orth;
		foreach my $id2 (keys(%$l2)){
			next if(!$f1_hits->{$id2});
			next if(!$f1_hits->{$id2}{$id1});
			next if(!$f2_hits->{$id1}{$id2});
			#print "\t[$id2]";
			#print "HERE!\t";
			$pairs{"$id1:$id2"}=$f1_hits->{$id2}{$id1} + $f2_hits->{$id1}{$id2};
			my @coverage = (($f1_lengths->{$id2}{$id1}{$id2}/$l2->{$id2}), ($f1_lengths->{$id2}{$id1}{$id1}/$l1->{$id1}), ($f2_lengths->{$id1}{$id2}{$id2}/$l2->{$id2}), ($f2_lengths->{$id1}{$id2}{$id1}/$l1->{$id1}));
			$coverage{"$id1:$id2"}=minA( @coverage  );
		}
	}
	
	
	my %seen;
	#Sort on sum of reciprocal phmmer bitscores:
	foreach my $p ( sort {$pairs{$b} <=> $pairs{$a}} keys %pairs){
		my ($id1, $id2) = split(/:/, $p);
		next if defined $seen{$id1};
		next if defined $seen{$id2};
		if ($pairs{$p} > 20 && $coverage{"$id1:$id2"}>0.75){
			$orths{$id1} = $id2;
			($seen{$id2},$seen{$id1})=(1,1);
			#print "$id1\t$id2\t$pairs{$p}\t" . $coverage{"$id1:$id2"} . "\n";
			#print "YES!\n";
		}
	}
	return \%orths;
}

sub find_seq_lengths {
	my($fasta) = @_;
	my ($name);
	my %l;
	open(STAT, "esl-seqstat -a $fasta |") or die "Unable to open pipe for [esl-seqstat -a $fasta].\n[$!]";
	while(<STAT>){
		if(/^\=\s+(\S+)\s+(\d+)/){
			$l{$1}=$2;
			#print "find_seq_lengths:\t[$1]\t[$2]\n";
		}
	}
	return(\%l);
}

sub parse_phmmer_tbl3 {
	my($tbl) = @_;
	my (%hits, %lengths);
	open TBL, "<", $tbl or die "Cannot open $tbl: $!";
	while(<TBL>){
		next if($_ =~ /^#/);
			chomp;
		my @splat = split /\s+/;
		$hits{$splat[0]}{$splat[3]} = 0 if (not defined($hits{$splat[0]}{$splat[3]}));
		$hits{$splat[0]}{$splat[3]} += $splat[7]; #really should check if hits overlap
		
		$lengths{$splat[0]}{$splat[3]}{$splat[0]} = 0 if (not defined($lengths{$splat[0]}{$splat[3]}{$splat[0]}));
		$lengths{$splat[0]}{$splat[3]}{$splat[3]} = 0 if (not defined($lengths{$splat[0]}{$splat[3]}{$splat[3]}));
		$lengths{$splat[0]}{$splat[3]}{$splat[3]} += ($splat[16] - $splat[15] + 1); #hmm coord -- query
		$lengths{$splat[0]}{$splat[3]}{$splat[0]} += ($splat[18] - $splat[17] + 1); #ali coord -- target
	}
	close TBL;
	return 0 if(! scalar(keys(%hits)));
	return (\%hits, \%lengths);
}

sub minA {
	my $min = $_[0];
	foreach my $a (@_){
		$min = min($min, $a) if isNumeric($a);
	}
	return $min;
}

sub isNumeric {
	my $num = shift;
	if ($num=~/^-?\d+\.?\d*$/) {
		return 1;
	}
	else {
		return 0;
	}
}

#max
sub max {
	return $_[0] if @_ == 1;
	$_[0] > $_[1] ? $_[0] : $_[1]
}

#min
sub min {
	return $_[0] if @_ == 1;
	$_[0] < $_[1] ? $_[0] : $_[1]
}
