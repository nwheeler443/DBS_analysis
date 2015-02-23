#!/usr/bin/perl
##!/software/perl-5.8.8/bin/perl

use warnings;
use strict;

use Bio::Perl;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::SeqFeature::Generic;
use Statistics::Distributions;
use Data::Dumper;

my $embl = shift @ARGV;
my $outfile = shift @ARGV;
my %gene_info;

my $cds_coordinates  = cds_locations($embl);
my $annotation_file =  Bio::SeqIO->new(-file => $embl, -format => 'EMBL') or die "Error: Couldnt open embl file: $!\n";
#my $annotation_file =  Bio::SeqIO->new(-file => $embl, -format => 'GENBANK') or die "Error: Couldnt open GENBANK file: $!\n";

open OUT, ">$outfile" or die "Couldn't open file for writing: $!\n";

while (my $sequence_annotation = $annotation_file->next_seq())      # work this out
{
    for my $feature ($sequence_annotation->get_SeqFeatures())       # get_SeqFeatures is part of a package?
    {
        next if !($feature->primary_tag eq 'CDS' || $feature->primary_tag eq 'polypeptide');
		#next if $feature->has_tag("pseudo");
		next if $feature->has_tag("pseudo") || $feature->has_tag("pseudogene");
        my $feature_id    = get_feature_id($feature);
        my $gene_name     = get_gene_name($feature);
        my $product_value = get_product_value($feature);
        my $start = $feature->start;
        my $end = $feature->end;
        my $sequence = '';
        
        @{$gene_info{$feature_id}} = ($gene_name, $product_value);
        if($feature->has_tag("translation")){
            my @value = $feature->get_tag_values("translation");
            $sequence = $value[0];
        }
		
		elsif($feature->strand == 1){
            $sequence = translate_as_string($feature->entire_seq->subseq($start, $end), -complete=>1,
            -codontable_id=>11);
        } else {
            $sequence = translate_as_string(revcom($feature->entire_seq->subseq($start, $end))->seq(), -complete=>1,
            -codontable_id=>11);
        }
        
        print OUT ">$feature_id\n";
        print OUT substr($sequence, 0, 80, '')."\n" while (length($sequence));
    }
}
close OUT;

sub cds_locations
{
    my($embl_file) = @_;
    my @cds_coordinates;
    
	my $annotation_file =  Bio::SeqIO->new(-file => $embl_file, -format => 'EMBL') or die "Error: Couldnt open the annotation file\n";
	#my $annotation_file =  Bio::SeqIO->new(-file => $embl_file, -format => 'GENBANK') or die "Error: Couldnt open the annotation file\n";
    while (my $sequence_annotation = $annotation_file->next_seq())
    {
        for my $feature ($sequence_annotation->get_SeqFeatures())
        {
            next if !($feature->primary_tag eq 'CDS');
            push(@cds_coordinates, [$feature->start,$feature->end]);
        }
    }
    return \@cds_coordinates;
}

sub get_feature_id
{
    my($feature) = @_;
    my $feature_id = int(rand(10000));
    my @junk;
    if($feature->has_tag('locus_tag'))
    {
        ($feature_id, @junk) = $feature->get_tag_values('locus_tag');
    }
    elsif($feature->has_tag('ID'))
    {
        ($feature_id, @junk) = $feature->get_tag_values('ID');
    }
    elsif($feature->has_tag('systematic_id'))
    {
        ($feature_id, @junk) = $feature->get_tag_values('systematic_id');
    }
    else
    {
        $feature_id = join("_",($feature->seq_id(), $feature->strand, $feature->start, $feature->end ));
    }
    $feature_id =~ s/^"|"$//g;
    return $feature_id ;
}

sub get_gene_name
{
    my($feature) = @_;
    my $gene_name;
    my @junk;
    if($feature->has_tag('gene'))
    {
        ($gene_name, @junk) = $feature->get_tag_values('gene');
    }
    else
    {
        $gene_name = get_feature_id($feature);
    }
    $gene_name =~ s/\W//g;
    return $gene_name;
}

sub get_product_value
{
    my($feature) = @_;
    my $product = "";
    my @junk;
    if($feature->has_tag('product'))
    {
        ($product, @junk) = $feature->get_tag_values('product');
    }
    
    return $product;
}

sub is_gene_within_cds
{
    my($cds_coordinates, $gene_feature) = @_;
    for my $current_coords(@{$cds_coordinates})
    {
        next if( $current_coords->[0] > $gene_feature->start);
        next if( $current_coords->[1] < $gene_feature->end);
        return 1;
    }
    
    return 0;
}