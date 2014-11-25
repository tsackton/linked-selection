#!/usr/bin/env perl

use strict;
use warnings;
use POSIX qw(ceil floor);

#Specify the path to the directory of files to process as the first argument
my $filepath = $ARGV[0];
die "Need a file path to gzipped VCF files as first argument" unless $filepath;

#Specify the window size in kb as the second argument
my $windowsize = $ARGV[1];

#Specify outfile prefix as the third argument
my $outfileprefix = $ARGV[2];

unless ($windowsize) {
	warn "No window size specificed as second argument. Using default 500kb window size.\n";
	$windowsize = 500000;
}

unless ($outfileprefix) {
	warn "Using default outfile prefix (window).\n";
	$outfileprefix="window";
}

my %windowstats;
my %avgpi;

my @files = <$filepath/*.gz>;

print "Processing files....\n";

foreach my $file (@files) {
	my $opencmd = "gzip -cd $file |";
	open IN, $opencmd or die;
	my $sp = $file;
	$sp =~ s/.+\/(\w+)\..+/$1/;
	print "Working on species $sp in file $file\n";	
	
	while (<IN>) {
		chomp;
		my ($chr, $pos, undef, undef, $pi) = split;
		#assumes VCF format is chr, position, chromosomes, minor allele freq, pi, separated by white space (space or tab)
		
		#does not error check to make sure that chromosome names do not contain spaces, so this needs to be confirmed prior to input
		
		##CLEAN UP CHROMOSOMES##	
		##the following code is specific to the set of organisms and is used to convert
		##the wide variety of chromosome naming schemes into simple integers
		##also removes sex chromosomes (X, Y, Z, W) and mt, if present
		##finally, this also removes minor scaffolds for species where the chromosomal assembly
		##includes both the major (chromosome-scale) scaffolds and the minor scaffolds	
		
		#fix Agam chromosomes
		if ($sp eq "Agam" and $chr eq "2L") {
			$chr = "2";
			$pos = $pos + 61545105;
		}
		elsif ($sp eq "Agam" and $chr eq "3L") {
			$chr = "3";
			$pos = $pos + 53200684;
		}
		elsif ($sp eq "Agam" and $chr eq "2R") {
			$chr = "2";
		}
		elsif ($sp eq "Agam" and $chr eq "3R") {
			$chr = "3";
		}
		
		#MAIN CLEANUP SUBROUTINE
		$chr = clean_chr($chr, $sp);
		next unless $chr;
	
		#END CHR CLEANING#
		
		my $window = ceil($pos / $windowsize); #Compute window index
		$windowstats{$sp}{$chr}{$window}{'pi'} += $pi; #Add pi to cumulative sum for window
		$windowstats{$sp}{$chr}{$window}{'sites'}++; #Count number of sites in window
		$avgpi{$sp}{'pi'} += $pi; #Average pi for species
		$avgpi{$sp}{'sites'}++; #Number of sites for species-wide average pi
	}
}

print "Calculating windows...\n";

my $kbsize = $windowsize / 1000;

my $outfile = "$outfileprefix.$kbsize.txt";

open OUT, ">$outfile" or die;
foreach my $sp (keys %windowstats) {
	foreach my $chr (keys %{$windowstats{$sp}}) {
		foreach my $window (sort {$a <=> $b} keys %{$windowstats{$sp}{$chr}}) {	
			my $pisum = $windowstats{$sp}{$chr}{$window}{'pi'}; #Total pi for window
			my $sites = $windowstats{$sp}{$chr}{$window}{'sites'}; #Total sites in window
			my $avgpi;
			if (!$sites) {
				#No sites found for window
				$sites = 0;
				$avgpi = "NA";
			}
			else {
				$avgpi = $pisum / $sites; #pi per site for window
			}
			my $windpos = ((($window-1)*$windowsize) + ($window*$windowsize))/2; #midpoint of window is used as the window position
			print OUT "$sp\t$chr\t$window\t$windpos\t$avgpi\t$sites\n"; #output is species, chromosome, window number, window position, pi per bp in the window, number of sites in the window
		}
	}
}

open AVG, ">$outfileprefix.average_pi.txt" or die; #overall average pi
foreach my $sp (keys %avgpi) {
	my $ave = $avgpi{$sp}{'pi'} / $avgpi{$sp}{'sites'};
	print AVG "$sp\t$ave\t$avgpi{$sp}{'sites'}\n";
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
	
	if ($sp eq "Bdistachyon" and $chr =~ /scaffold/) {
		$chr = "";
		return ($chr);
	}

	if ($sp eq "Bmandarina" and $chr eq "chr1") {
		$chr = "";
		return($chr);
	}

	if ($chr eq "X") {
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
		$chr = "" if $chr eq "4";
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

	return $chr;		
}
