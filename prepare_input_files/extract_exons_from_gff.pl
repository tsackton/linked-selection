#!/usr/bin/env perl

use strict;
use warnings;

#extracts coding exons (CDS) from a GFF/GTF file and makes a BED file
#requires GTF/GFF format and defined CDS features
#processes all files in a directory that have an extension that starts with g
#output is a .exon file with the same prefix as the GTF/GFF file

my @files = <*.g*>;

foreach my $file (@files) {
	my $sp = $file;
	$sp =~ s/\.\w+$//;
	open IN, $file or die;
	open OUT, ">$sp.exon" or die;
	while (<IN>) {
		chomp;
		my @fields = split("\t", $_);
		next unless defined($fields[2]);
		next unless $fields[0];
		next unless $fields[2] eq "CDS";
		$fields[3] = $fields[3] - 1;
		my $start;
		my $end;
		if ($fields[3] > $fields[4]) {
			$start = $fields[4];
			$end = $fields[3];
		}
		else {
			$start = $fields[3];
			$end = $fields[4];
		}
		print OUT "$fields[0]\t$start\t$end\n";
	}
}