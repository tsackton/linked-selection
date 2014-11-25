#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(ceil floor);

my $filepath = $ARGV[0];
die "Need a file path to gzipped VCF files as first argument" unless $filepath;

#Specify the window size in kb as the second argument
my $windowsize = $ARGV[1];

unless ($windowsize) {
	warn "No window size specificed as second argument. Using default 500kb window size.\n";
	$windowsize = 500000;
}

my @files = <$filepath/*.exon>;
my %gd; #hash to store gene density in windows
my %total; #hash to store total coding bases per species

print "Processing files....\n";

foreach my $file (@files) {
	my $sp = $file;
	$sp =~ s/.+\/(\w+)\..+/$1/;
	print "Working on species $sp in file $file\n";	

	my $newfile = "${sp}.fix";
	my $fixoverlap = "sortBed -i $file | mergeBed -i - > ${sp}.fix"; 
	system($fixoverlap) unless (-f $newfile);
	open IN, "$newfile" or die;
	while (<IN>) {
		chomp;
		my ($chr, $start, $end) = split;
		$start = $start+1;
		
		##CLEAN UP CHROMOSOMES##	
		##the following code is specific to the set of organisms and is used to convert
		##the wide variety of chromosome naming schemes into simple integers
		##also removes sex chromosomes (X, Y, Z, W) and mt, if present
		##finally, this also removes minor scaffolds for species where the chromosomal assembly
		##includes both the major (chromosome-scale) scaffolds and the minor scaffolds	
		
		#fix Agam chromosomes
		if ($sp eq "Agam" and $chr eq "2L") {
			$chr = "2";
			$start = $start + 61545105;
			$end = $end + 61545105;
		}
		elsif ($sp eq "Agam" and $chr eq "3L") {
			$chr = "3";
			$start = $start + 53200684;
			$end = $end + 53200684;
		}
		elsif ($sp eq "Agam" and $chr eq "2R") {
			$chr = "2";
		}
		elsif ($sp eq "Agam" and $chr eq "3R") {
			$chr = "3";
		}
		
		$total{$sp} += $end-$start;
		
		#MAIN CLEANUP SUBROUTINE
		$chr = clean_chr($chr, $sp);
		next unless $chr;
	
		#END CHR CLEANING#
		
		for my $pos ($start..$end) {
			my $window = ceil($pos / $windowsize);
			$gd{$sp}{$chr}{$window}++;
		}
	}
	close IN;
	system("rm $newfile");
}

print "Calculating windows...\n";

my $kbsize = $windowsize / 1000;

open OUT, ">gd.$kbsize.txt" or die;
foreach my $sp (keys %gd) {
	foreach my $chr (keys %{$gd{$sp}}) {
		foreach my $window (sort {$a <=> $b} keys %{$gd{$sp}{$chr}}) {	
			my $exonsum = $gd{$sp}{$chr}{$window} || 0;
			my $length = $windowsize;
			my $windpos = ((($window-1)*$windowsize) + ($window*$windowsize))/2;
			my $exonfrac = $exonsum / $length;
			print OUT "$sp\t$chr\t$window\t$windpos\t$exonsum\t$exonfrac\n";
		}
	}
}

open TOTAL, ">total.exons.txt" or die;
foreach my $sp (keys %total) {
	my $exonbp = $total{$sp};
	print TOTAL "$sp\t$exonbp\n";
}
sub clean_chr {
	my %roman_convert = qw(I 1 II 2 III 3 IV 4 V 5 VI 6 VII 7 VIII 8 IX 9 X 10 XI 11 XII 12 XIII 13 XIV 14 XV 15 XVI 16 XVII 17 XVIII 18 XIX 19 XX 20 XXI 21 XXII 22 XXIII 23 XXIV 24 XXV 25);
	my ($chr, $sp) = @_;
	
	$chr =~ s/^CHROMOSOME_//;
	$chr =~ s/^chromosome_//;
	$chr =~ s/^chr//;
	$chr =~ s/^Chr//;
	$chr =~ s/lg//;
	$chr =~ s/Bd//;

	if (($sp eq "Celegans" or $sp eq "Cbriggsae") and ($chr =~ /X/)) {
		$chr = "";
		return($chr);
	}
	
	if (exists($roman_convert{$chr}) and ($sp eq "Gaculeatus" or $sp eq "Celegans" or $sp eq "Cbriggsae")) {
		$chr = $roman_convert{$chr};
	}
	
	if ($chr eq "19" and $sp eq "Gaculeatus") {
		$chr="";
		return($chr);
	}
	
	if ($chr eq "UNKN") {
		$chr = "";
		return ($chr);
	}
	
	if ($chr eq "M") {
		$chr = "";
		return($chr);
	}
	
	if ($chr eq "C") {
		$chr = "";
		return ($chr);
	}
	
	if ($sp eq "Bdistachyon" and $chr =~ /scaffold/) {
		$chr = "";
		return ($chr);
	}

	if ($sp eq "Bmandarina" and $chr eq "1") {
		$chr = "";
		return($chr);
	}

	if ($chr eq "X") {
		$chr = "";
		return($chr);
	}
	
	if ($chr eq "0") {
		$chr = "";
		return($chr);
	}
	
	if ($chr eq "MtDNA") {
		$chr = "";
		return($chr);
	}
	
	if ($chr =~ /random/) {
		$chr = "";
		return($chr);
	}
	
	if ($chr eq "MT") {
		$chr = "";
		return($chr);
	}
	
	if ($chr eq "Z") {
		$chr = "";
		return($chr);
	}
	
	if ($sp eq "Clupus" and $chr =~ /[A-Z]+\d+\.\d+/) {
		$chr = "";
		return($chr);
	}
	
	if ($sp eq "Cret") {
		$chr =~ s/^scaffold_//;
		$chr = "" if  $chr > 9;
		return ($chr);
	}
	
	if ($sp eq "Crubella") {
		$chr =~ s/^scaffold_//;
		$chr = "" if  $chr > 8;
		return ($chr);
	}
		
	if ($sp eq "Csativus") {
		$chr = "" if $chr =~ /Scaffold/;
		return($chr);
	}
	
	if ($sp eq "Dmelanogaster") {
		$chr = "" unless ($chr eq "3" or $chr eq "5" or $chr eq "7" or $chr eq "8");
		return ($chr);
	}
		
	if ($sp eq "Dpseudoobscura") {
		$chr = "" unless $chr eq "2";
		return ($chr);
	}
	
	if ($sp eq "Drerio") {
		$chr = "" if $chr =~ /^Zv9_/;
		return($chr);
	}
	
	if ($chr eq "W") {
		$chr = "";
		return($chr);
	}
	
	if ($chr eq "Y") {
		$chr = "";
		return($chr);
	}
	
	if ($sp eq "Gmax") {
		$chr =~ s/Gm//;
		$chr = "" if $chr =~ /scaffold/;
		return($chr);
	}
	
	if ($sp eq "Gaculeatus") {
		$chr = "" if $chr =~ /scaffold/;
		return($chr);
	}
	
	if ($sp eq "Graimondii") {
		$chr =~ s/^0//;
		$chr = "" if $chr =~ /scaffold/;
	}
	
	if ($sp eq "Loculatus") {
		$chr =~ s/LG//;
		$chr = "" if $chr =~ /[A-Z]+\d+\.\d+/;
		return($chr);
	}
	
	if ($chr =~ /PATCH/) {
		$chr = "";
		return($chr);
	}
	
	if ($sp eq "Mmul" and $chr =~ /^NW_/) {
		$chr = "";
		return($chr);
	}
	
	if ($sp eq "Olatipes") {
		$chr = "" if $chr =~ /scaffold/;
		$chr = "" if $chr =~ /ultracontig/;
		return($chr);
	}
	
	if ($sp eq "Mmcastaneus" ) {
		$chr = "" if $chr =~ /[A-Z]+\d+\.\d+/;
		return($chr);
	}
	
	if ($chr eq "Sy") {
		$chr = "";
		return($chr);
	}
	
	if ($sp eq "Pdav" and $chr > 8) {
		$chr = "";
		return($chr);
	}
	
	if ($sp eq "Ptrichocarpa" and $chr =~ /scaffold/) {
		$chr = "";
		return ($chr);
	}
	
	if ($chr eq "Un") {
		$chr = "";
		return($chr);
	}
	
	if ($sp eq "Sbicolor" and $chr =~ /super/) {
		$chr = "";
		return($chr);		
	}
	
	if ($sp eq "Sitalica") {
		$chr =~ s/scaffold_//;
		$chr = "" if $chr > 9;
		return($chr);
	}
	
	if ($chr eq "UNKNOWN") {
		$chr = "";
		return($chr);
	}
	
	if ($chr eq "Pt") {
		$chr = "";
		return($chr);
	}
	if ($sp eq "Sscrofa" and $chr =~ /[A-Z]+\d+/) {
		$chr = "";
		return($chr);
	}
	
	if ($sp eq "Csemilaevis") {
		my %csemkey = qw(NC_024307.1 1 NC_024308.1 2 NC_024309.1 3 NC_024310.1 4 NC_024311.1 5 NC_024312.1 6 NC_024313.1 7 NC_024314.1 8 NC_024315.1 9 NC_024316.1 10 NC_024317.1 11 NC_024318.1 12 NC_024319.1 13 NC_024320.1 14 NC_024321.1 15 NC_024322.1 16 NC_024323.1 17 NC_024324.1 18 NC_024325.1 19 NC_024326.1 20);
		if (exists($csemkey{$chr})) {
			$chr = $csemkey{$chr};
		}
		else {
			$chr = "";
		}
		return ($chr);
	}

	if ($chr eq "un") {
		$chr = "";
		return($chr);
	}
	
	if ($chr eq "Mt") {
		$chr = "";
		return($chr);
	}

	return $chr;		
}
