#!/usr/bin/perl 


use warnings;
use strict;

# FT                   /gene="thrA"
# FT                   /locus_tag="KPN_00002"
# FT                   /product="bifunctional aspartokinase I/homeserine
# FT                   dehydrogenase I"

# FT   CDS             complement(14251..14997)
# FT                   /codon_start=1
# FT                   /transl_table=11
# FT                   /locus_tag="KPK_0012"
# FT                   /product="transcriptional regulator, GntR family"



my $fileName = shift; 

my ($printMe,$tag,$name,$product)=(0,"","","");
open(F, "< $fileName");
while (my $f = <F>){
    chomp($f);
    if ($f=~/^FT\s+\/gene="(\S+)"/){
	$name=$1;	
	$printMe=0;
    }
    elsif($f=~/^FT\s+\/locus_tag="(\S+)"/){
	$tag=$1; 	
    }
    elsif($f=~/^FT\s+\/product="(.*)"/){
	$product=$1; 
	$printMe=1; 
    }
    elsif($f=~/^FT\s+\/product="(.*)/){
	$product=$1; 
	
	while(!$printMe){
	    $f = <F>; 
	    chomp($f); 
	    if($f=~/^FT\s+(.*)"/){
		$product.=$1; 
		$printMe=1; 
	    }
	    else{
		$product.=$1; 
	    }
	}
    
    }
    
    if($printMe){
	print "$tag\t$name\t$product\n"; 
	($printMe,$tag,$name,$product)=(0,"","","");
    }
    
}
close(F);


exit(0); 



