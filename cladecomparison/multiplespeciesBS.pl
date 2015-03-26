#!/usr/bin/perl


# have decided to assign one score per domain per protein per isolate - if there are multiple copies I will only look at one

use warnings;
use strict;
use Data::Dumper;

# hash of all genes from each GS and the associated domains
my $ref1 = "REF1";						# reference form first group
my $ref2 = "REF2";						# reference from second group

my %gs1arch;
my %gs2arch;

open DBS, "$ref1-$ref2/results.dbs";
while (<DBS>) {
	if ($_ =~ /($ref1\S+)\s+(\S+)\s+(PF\S+)\s+\d/) {
		push @{$gs1arch{$1}}, $3;
		push @{$gs2arch{$2}}, $3;
	}
}

close DBS;

# retrieve scores for the representatives - made a hash (gene names) of hashes (domain names)

my %GS1scores;

open IN, "genomes/$ref1-pfam_hmmscan1.tbl";
while (<IN>) {
	next if ($_ =~ /\#/);
	chomp;
	my @split = split(/\s+/, $_);
	if ($split[1] ~~ @{$gs1arch{$split[3]}}) {
		@{$GS1scores{$split[3]}{$split[1]}}[0] = $split[7];		# some of these proteins will have multiple instances of the same domain
	}
}
close IN;

# for each isolate, find orthologs and pair them with domain architectures, then retrieve these from the hmmscan file

my @GS1 = ("CCUG19995", "UNSWCD", "UNSW3", "UNSWCS", "UNSW1", "Lasto127", "13826", "ATCC51561");
my @GS2 = ("Lasto64", "ATCC51562", "Lasto220", "Lasto393", "Lasto28", "ATCC33237", "Lasto205");

# make hash of orthologs, then for each entry in hmmscan, match gene to orth then search for orth in gs1arch
foreach my $number (0..$#GS1) {
	my $name = $GS1[$number];
	$number++;
	my %orthlist;
	open IN, "orthlists/$ref1-$name";
	while (<IN>) {
		if ($_ =~ /(\S+)\s+(\S+)/) {
			$orthlist{$2} = $1;
		}
	}
	close IN;
	open IN, "genomes/$name-pfam_hmmscan1.tbl";
	while (<IN>) {
		next if ($_ =~ /\#/);
		chomp;
		my @split = split(/\s+/, $_);
		if ($split[1] ~~ @{$gs1arch{$orthlist{$split[3]}}}) {		# if the domain is associated with the ortholog to that gene
			@{$GS1scores{$orthlist{$split[3]}}{$split[1]}}[$number] = $split[7];
		}
	}
	close IN;
}

#print Dumper (\%GS1scores);

my %GS2scores;

open IN, "genomes/$ref2-pfam_hmmscan1.tbl";
while (<IN>) {
	next if ($_ =~ /\#/);
		chomp;
	my @split = split(/\s+/, $_);
	if ($split[1] ~~ @{$gs2arch{$split[3]}}) {
		@{$GS2scores{$split[3]}{$split[1]}}[0] = $split[7];
	}
}
close IN;

# make hash of orthologs, then for each entry in hmmscan, match gene to orth then search for orth in gs1arch
foreach my $number (0..$#GS2) {
	my $name = $GS2[$number];
	$number++;
	my %orthlist;
	open IN, "orthlists/$ref2-$name";
	while (<IN>) {
		if ($_ =~ /(\S+)\s+(\S+)/) {
			$orthlist{$2} = $1;
		}
	}
	close IN;
	open IN, "genomes/$name-pfam_hmmscan1.tbl";
	while (<IN>) {
		next if ($_ =~ /\#/);
			chomp;
		my @split = split(/\s+/, $_);
		if ($split[1] ~~ @{$gs2arch{$orthlist{$split[3]}}}) {		# if the domain is associated with the ortholog to that gene
			@{$GS2scores{$orthlist{$split[3]}}{$split[1]}}[$number] = $split[7];
		}
	}
	close IN;
}

#print Dumper (\%GS2scores);

# print results to file!!

open OUT, "> GS1scores.tsv";
foreach my $gene (keys(%GS1scores)) {
	foreach my $domain (keys(%{$GS1scores{$gene}})) {
		print OUT "$gene\t$domain\t";
		foreach (@{$GS1scores{$gene}{$domain}}) {
			print OUT "$_\t";
		}
		print OUT "\n";
	}
}
close OUT;

open OUT, "> GS2scores.tsv";
foreach my $gene (keys(%GS2scores)) {
	foreach my $domain (keys(%{$GS2scores{$gene}})) {
		print OUT "$gene\t$domain\t";
		foreach (@{$GS2scores{$gene}{$domain}}) {
			print OUT "$_\t";
		}
		print OUT "\n";
	}
}

close OUT














