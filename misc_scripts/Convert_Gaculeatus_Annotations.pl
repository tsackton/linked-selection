use strict ;
use warnings ; 

=begin
 
 This script will read in the original gtf file for Gaculateaus and convert the coordinates to the improved chromosome scale assembly as described in :
 Recombination in the threespine stickleback genome – patterns and consequences. Molecular Ecology 22(11): 3014–3027.
 
=cut

#### read in gtf file line by line and convert on teh fly
while (<STDIN>) { 
	chomp $_ ; 
	my @split = split ( /\t/, $_ ) ; 
	my ( $seq, $coord1, $strand ) = convert_coords( $split[0], $split[3], $split[6] ) ;
	my ( $xxx, $coord2, $yyy ) = convert_coords( $split[0], $split[4], $split[6] ) ; 
	print $seq, "\t", $split[1], "\t", $split[2], "\t", $coord1, "\t", $coord2, "\t", $split[5], "\t", $strand, "\t", $split[7], "\t", $split[8], "\n" ;
}

### this portion will convert coordinates from the original assembly to the new coordiantes
sub convert_coords {

	my ($seqname, $coord, $strand) = @_ ;

	#first deal with flips	
	if ($seqname eq "chrI" && ($coord >= 19922392 && $coord <= 27290228) ) {
		$coord = 27290228 - ($coord - 19922392);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "chrI" && ($coord >= 27291229 && $coord <= 28185914) ) {
		$coord = 28185914 - ($coord - 27291229);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "chrIV" && ($coord >= 17820057 && $coord <= 28355687)) {
		$coord = 28355687 - ($coord - 17820057);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "chrV" && ($coord >= 1 && $coord <= 7877689)) {
		$coord = 7877689 - ($coord - 1);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "chrV" && ($coord >= 10998358 && $coord <= 11741737)) {
		$coord = 11741737 - ($coord - 10998358);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "chrVI" && ($coord >= 1 && $coord <= 4056779)) {
		$coord = 4056779 - ($coord - 1);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "chrIX" && ($coord >= 4309200 && $coord <= 17894190)) {
		$coord = 17894190 - ($coord - 4309200);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}	
	elsif ($seqname eq "chrXII" && ($coord >= 2209412 && $coord <= 16547614)) {
		$coord = 16547614 - ($coord - 2209412);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "chrXIII" && ($coord >= 18482172 && $coord <= 19650590)) {
		$coord = 19650590 - ($coord - 18482172);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "chrXIV" && ($coord >= 13083936 && $coord <= 15246461)) {
		$coord = 15246461 - ($coord - 13083936);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "chrXVI" && ($coord >= 16169748 && $coord <= 18115788)) {
		$coord = 18115788 - ($coord - 16169748);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}		
	elsif ($seqname eq "chrXIX" && ($coord >= 3824254 && $coord <= 20240660)) {
		$coord = 20240660 - ($coord - 3824254);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}	
	elsif ($seqname eq "chrXX" && ($coord >= 834142 && $coord <= 17950836)) {
		$coord = 17950836 - ($coord - 834142);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}

	#then add on scaffolds that are at the right periphery
	if ($seqname eq "scaffold_90") {
		$seqname = "chrIV";
		$coord += (32632949 - 1) + (33024843 - 32632949 + 1);
	}
	elsif ($seqname eq "scaffold_121") {
		$seqname = "chrIV";
		$coord = 33024843 - ($coord - 1);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "scaffold_112") {
		$seqname = "chrV";
		$coord = 12609052 - ($coord - 1);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "scaffold_156") {
		$seqname = "chrVI";
		$coord = 12251398 + $coord - 1;
	}
	elsif ($seqname eq "scaffold_69") {
		$seqname = "chrXVI";
		$coord += (18115789 - 1) + (18459807 - 18115789 + 1);
	}
	elsif ($seqname eq "scaffold_115") {
		$seqname = "chrXVI";
		$coord = 18459807 - ($coord - 1);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "scaffold_67") {
		$seqname = "chrXXI";
		$coord = 11717488 + $coord - 1;
	}

	#then deal with shifting because of insertions
	if ($seqname eq "chrI" && $coord >= 27290229) {
		$coord += (27972651 - 27290229 + 1);
	}
	elsif ($seqname eq "chrV" && $coord >= 1) {
		$coord += ((1738956 - 1 + 1) + (1212556 - 1 + 1));
	}
	elsif ($seqname eq "chrVI" && $coord >= 1) {
		$coord += (1001911 - 1 + 1);
	}
	elsif ($seqname eq "chrVII" && $coord >= 10410930) {
		$coord += (12159109 - 10410930 + 1);
	}
	elsif ($seqname eq "chrX" && $coord >= 1) {
		$coord += (1110135 - 1 + 1);
	}
	elsif ($seqname eq "chrXI" && $coord >= 15152914) {
		$coord += (15344290 - 15152914 + 1);
	}
	elsif ($seqname eq "chrXII" && $coord >= 1247383) {
		$coord += (2116554 - 1247383 + 1);
	}
	elsif ($seqname eq "chrXVII" && $coord >= 6449158) {
		$coord += (11534418 - 6449158 + 1);
	}
	elsif ($seqname eq "chrXX" && $coord >= 1) {
		$coord += (304788 - 1 + 1);
	}
	elsif ($seqname eq "chrXXI" && $coord >= 1) {
		$coord += (2648413 - 1 + 1);
	}
	
	#finally deal with converting remaining scaffolds to correct chromosome locations
	
	if ($seqname eq "scaffold_74") {
		$seqname = "chrI";
		$coord += 27290229 - 1;
	}
	elsif ($seqname eq "scaffold_54") {
		$seqname = "chrV";
	}
	elsif ($seqname eq "scaffold_48") {
		$seqname = "chrV";
		$coord += 1212556;
	}
	elsif ($seqname eq "scaffold_61") {
		$seqname = "chrVI";
	}
	elsif ($seqname eq "scaffold_47") {
		$seqname = "chrVII";
		$coord += 10410930 - 1;
	}
	elsif ($seqname eq "scaffold_58") {
		$seqname = "chrX";
	}
	elsif ($seqname eq "scaffold_151") {
		$seqname = "chrXI";
		$coord += 15152914 - 1;
	}
	elsif ($seqname eq "scaffold_68") {
		$seqname = "chrXII";
		$coord += 1247383 - 1;
	}
	elsif ($seqname eq "scaffold_27") {
		$seqname = "chrXVII";
		$coord = 11534418 - ($coord - 1);
		$strand = $strand eq "+" ? $strand = "-" : $strand = "+";
	}
	elsif ($seqname eq "scaffold_122") {
		$seqname = "chrXX";
	}
	elsif ($seqname eq "scaffold_37") {
		$seqname = "chrXXI";
	}
	
	return ($seqname, $coord, $strand);
}
	
