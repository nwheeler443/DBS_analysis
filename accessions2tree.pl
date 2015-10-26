#!/usr/bin/perl 

use warnings;
use strict; 

my @accessions=@ARGV; 

print "[@accessions]\n"; 

my (%acc2name,%name2acc,%nameCnt); 
foreach my $acc (@accessions){
    
    my $sum = `grep ^$acc summary.txt`; 
    next if (not defined($sum)); 
    my @sum=split(/\t/, $sum); 
    my $shortname = species2shortspecies($sum[6]);
    
    if (not defined($nameCnt{$shortname})){
	$nameCnt{$shortname}=1; 
	$name2acc{$shortname}=$acc;
	$acc2name{$acc}=$shortname;
    }
    else{
	$nameCnt{$shortname}++;
	$shortname.=$nameCnt{$shortname};
	$name2acc{$shortname}=$acc;
	$acc2name{$acc}=$shortname;	
    }
    print "$acc\t$shortname\t[$sum[6]]\n";
}

my $egrepStr = join('|', @accessions); 

open(IN, "egrep \'$egrepStr\' bacteria_RF00177.stk | grep -v \'^\#' | "); 
#esl-reformat pfam bacteria_RF00177.stk | egrep 'U000|AL12345|AP006716|AJ938182' | perl -lane 'if(/^#=GS (\S+)\.\d+\/\d+\-\d+/){$sum=`grep $1 /data/genomes/bacteria/summary.txt`; @sum=split(/\t/, $sum); print $sum[6]}'
open(OUT, "> tmpOut.txt");
print OUT "# STOCKHOLM 1.0\n\n";
while(my $in=<IN>){
    
    if($in=~/^(\S+)\.\d+\/\d+\-\d+\s+(\S+)/){
	my ($acc,$seq) = ($1,$2);
	next if (not defined($acc2name{$acc}));
	print OUT "$acc2name{$acc}    $seq\n";
    }
    
    print $in;
    
}
print OUT "\/\/\n";
close(OUT);

system("esl-reformat phylip  tmpOut.txt > infile && rm -f outfile outtree && echo \"Y\" | ./dnaml "); 



exit(0);
######################################################################

######################################################################
#species2shortspecies: Given a species string eg. "Homo sapiens
#                      (human)" generate a nicely formated short name
#                      with no whitespace eg. "H.sapiens".
sub species2shortspecies {
    my $species = shift;
    my $shortSpecies;
    
    if ($species=~/(.*)\s+sp\./){
	$shortSpecies = $1;
    }
    elsif ($species=~/metagenome/i or $species=~/uncultured/i){
	$species=~s/metagenome/metag\./gi;
	$species=~s/uncultured/uncult\./gi;
	my @w = split(/\s+/,$species);
	if(scalar(@w)>2){
	    foreach my $w (@w){
		$shortSpecies .= substr($w, 0, 5) . '.';
	    }
	}
	else {
	    $shortSpecies = $species;
	    $shortSpecies =~ s/\s+/_/g;
	}
    }#lots of conditions here. Need else you get some ridiculous species names.
    #elsif($species=~/^(\S+)\s+(\S{4,})/ && $species!~/[\/\-\_0-9]/ && $species!~/^[a-z]/ && $species!~/\svirus$/ && $species!~/\svirus\s/ && $species!~/^Plasmid\s/i && $species!~/\splasmid\s/i){
    elsif($species=~/^(\S+)\s+(\S{4,})/ && $species!~/^[a-z]/ && $species!~/\svirus$/ && $species!~/\svirus\s/ && $species!~/^Plasmid\s/i && $species!~/\splasmid\s/i){
	$shortSpecies = substr($1,0,1) . "." . $2; 
    }
    else {
	$shortSpecies = $species;
    }
    
    $shortSpecies =~ s/\s+/_/g;
    $shortSpecies =~ s/[\'\(\)\:\/]//g;
    $shortSpecies = substr($shortSpecies,0,9) if (length($shortSpecies) > 9);
    
#   H.P 
    return $shortSpecies;
}







