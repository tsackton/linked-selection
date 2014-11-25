#!/usr/bin/perl
use strict;
use warnings;

### This script converts positions from scaffold to chromosomes for Heliconius melpomene

### these hashes store relative positions for converting from each scaffold to chromosome
my %map_chr ;
my %map_fun ; 
my %map_coord ;

### this file is produced by script XXX
open IN, "<hmel_coord_conversion_table.txt" ;
while (<IN>) { 
	chomp $_ ; 
	my @split = split ( /\t/, $_ ) ; 
    $map_chr{"$split[0]"} = $split[1] ;     ### this hash stores the scaffold to chromosome relationship
    $map_fun{"$split[0]"} = 1 ;             ### this one stores whether coordinates are added or subtracted
	if ( $split[2] eq "subtract" ) { 
		$map_fun{"$split[0]"} = -1 ; 
	}
    $map_coord{"$split[0]"} = $split[3] ;   #### this stores the absolute value of coordinate differences
}
close IN ; 

### convert coordinates for snp files
### format of these files is
### scaffold /// position /// number of chromosomes represented /// MAF /// PI
while (<STDIN>) { 

    chomp $_ ;
	my @split = split ( /\t/, $_ ) ;
    if ( exists( $map_chr{$split[0]} ) ) {  ####
		my $coord_start = $map_coord{$split[0]} + $map_fun{$split[0]}*$split[3] ;
		my $coord_end = $map_coord{$split[0]} + $map_fun{$split[0]}*$split[4] ;
		print $map_chr{$split[0]}, "\t", $split[1], "\t", $split[2], "\t", $split[3], "\t", $split[4], "\n" ;
	}
}
