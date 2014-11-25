#!/usr/bin/env perl

use warnings ; 
use strict ; 
use Parallel::ForkManager;

#move files
system("cp ../rec_rate_est/*forgk.out.gz .");
system("cp ../rec_rate_est/mut_ests.txt . ");

#unzip
system("gunzip *.gz");

#make output dir
my $outdir = $ARGV[0];
system("mkdir -p $outdir");

my %mut; 
my @species;
open IN, "<mut_ests.txt" ;
while (<IN>) { 
	chomp $_ ; 
	next if /^species/;
	my @split = split ( /\t/, $_ ) ; 
	$mut{$split[0]}{'min'} = $split[6];
	$mut{$split[0]}{'const'} = $split[7];
	$mut{$split[0]}{'max'} = $split[8];
	
	push @species, $split[0]; 
}
close IN ; 

my @files = <*forgk.out>; 

my $manager = new Parallel::ForkManager( 60 );
foreach my $file ( @files ) { 
	foreach my $species (@species) {
		foreach my $mutkey (keys %{$mut{$species}}) {
			$manager->start and next;
			my $wind = $file;
			$wind =~ s/_forgk\.out//;
			my $outfile = "${species}.${wind}.${mutkey}.gk" ; 
			my $u = $mut{$species}{$mutkey};
			print "./compute_gk -s $species -u $u < $file > $outdir$outfile\n" ;
			system("./compute_gk -s $species -u $u < $file > $outdir$outfile") ; 		
			$manager->finish;
		}
	}
}
$manager->wait_all_children;

#clean up files
system("rm *forgk.out");