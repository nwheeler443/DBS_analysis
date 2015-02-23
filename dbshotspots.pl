#!/usr/bin/perl

use warnings;
use strict;
use Data::Dumper;

# input: ~/scripts/dbshotpots.pl commensal pathogen high/low

# read in dbs file and pick out a hash of domains with high dbs

my $comm = shift @ARGV;
my $path = shift @ARGV;
my $mode = shift @ARGV;

my %domainsofinterest;
open IN, "< $comm-$path/results.dbs" or die "couldn't open file";
while (<IN>) {
	# SEN0058	SG0061	PF03390.10	22	441	523.4	22	378	408.0	115.4	48.2413909578548	0
	if ($_ =~ /(\w+)\t(\w+)\t(PF\w+\.\d+)\t(\d+)\t(\d+)\t.+\t(\d+)\t(\d+)\t.+\t(\d+)\t/){
		#print $8;
		if ($mode eq "high") {
			if ($8 > 4) {
				#print $_;
				@{$domainsofinterest{$3}} = ($1, $4, $5, $2, $6, $7);
			}
		}
		if ($mode eq "low") {
			if ($8 < 4 && $8 > (-4
				)) {
				#print $_;
				@{$domainsofinterest{$3}} = ($1, $4, $5, $2, $6, $7);
			}
		}
	}
	
}
close IN;

#print Dumper (\%domainsofinterest);

# retrieve sequence data and align it

system("esl-sfetch --index genomes/$comm.fasta");
system("esl-sfetch --index genomes/$path.fasta");

# make a hash of alignments, with the domain name as the key

my %alignmentstring;

foreach my $domain (keys(%domainsofinterest)) {
	#print "$domainsofinterest{$domain}[0]\t$domainsofinterest{$domain}[3]\n";
	my $seq1 = `esl-sfetch -c $domainsofinterest{$domain}[1]..$domainsofinterest{$domain}[2] genomes/$comm.fasta $domainsofinterest{$domain}[0]`;
	my $seq2 = `esl-sfetch -c $domainsofinterest{$domain}[4]..$domainsofinterest{$domain}[5] genomes/$path.fasta $domainsofinterest{$domain}[3]`;
	#print "$seq1\t$seq2\n";
	open OUT, "> /tmp/sequences.txt";
	print OUT "$seq1$seq2";
	close OUT;
	$alignmentstring{$domain} = `mafft /tmp/sequences.txt`;
	system("rm /tmp/sequences.txt");
}

#print Dumper (\%alignments);

my %alignments;

foreach my $domain (keys (%alignmentstring)) {
	my @split = split (">\.+\n", $alignmentstring{$domain});
	@{$alignments{$domain}} = @split;
	
#	print Dumper (\%alignments);
}

my $mismatchcount = 0;
my %hits;

foreach my $domain (keys (%alignments)) {
	my @seq1 = split ("", $alignments{$domain}[1]);
	my @seq2 = split ("", $alignments{$domain}[2]);
	foreach my $start (0..$#seq1-25) {
		foreach my $position ($start..$start+25) {
			if ($seq1[$position] ne $seq2[$position]) {
				$mismatchcount ++;
			}
		}
		if ($mismatchcount > 5) {
			$hits{$domain}=1;
		}
	$mismatchcount = 0;
	}
}

foreach my $hit (keys(%hits)) {
	print ">", $hit, "_", $domainsofinterest{$hit}[0], "_species1\n", $alignments{$hit}[1], ">", $hit, "_species2\n", $alignments{$hit}[2], "\n";
}

