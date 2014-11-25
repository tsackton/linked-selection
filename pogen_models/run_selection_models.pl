#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;

#Perl script to automate analysis of selection models in R

my $indir = shift;
my @files = <$indir/*.gk>;

my $manager = new Parallel::ForkManager( 16 );
foreach my $file (@files) {
	$manager->start and next;
	$file =~ s/^.*?([A-Za-z._0-9]+)$/$1/;
	print "Processing $file...\n";
	my ($sp, $window, $rec, $umut, undef) = split(/[\._]/, $file);
	open my $params, ">$file.params" or die;
	$window =~ s/wind//;
	print $params "file = \"$indir\/$file\"\nspec = \"$sp\"\nwindsize = \"$window\"\nrecmethod = \"$rec\"\nuversion = \"$umut\"\n\n";
	print "R --quiet --no-save < $file.run\n";
	system("cat $file.params run_bgs.R > $file.run");
	system("R --quiet --no-save < $file.run");	
	$manager->finish;
}
$manager->wait_all_children
