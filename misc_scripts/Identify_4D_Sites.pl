#!/usr/bin/perl
use warnings ;
use strict ;

=begin

 This script takes a gff3 (or gtf or gff) file and a fasta file as input. 
 It then identifed all four-fold degenerate positions and retains those.
 
=cut

### file information is passed via command line
my $gff = $ARGV[0] ; 
my $fa = $ARGV[1] ; 

### read in set of 4d degenerate codons
### this file can take any format as long as all it contains are \t and 4-fold degenerate codons in the second-nth columns on each line
my %codon ; 
open CODON, "<codon_table.txt" ;
while(<CODON>) { 
	chomp $_ ;
	my @split = split ( /\t/, $_ ) ; 
	foreach ( 1..$#split ) { 
		$codon{$split[$_]} ++ ;
	}
}
close CODON ; 

### gene information is extracted from standard gff3 files.
### this also shoudl work for gtf and gff files
my %gene ; 
open GFF, "<$ARGV[0]" ;
while (<GFF>) { 
	chomp $_ ; 
	my @split = split ( /\t/, $_ ) ; 	
	### parse information for mRNA molecule
	if ( $split[2] eq "CDS" ) { 
		my $id = $split[0]."_".$split[3]."_".$split[4]."_".$split[7] ;
		$gene{$id}{"strand"} = $split[6] ;
		$gene{$id}{"chrom"} = $split[0] ; 
		$gene{$id}{"start"} = $split[3] - 1 ;
		$gene{$id}{"stop"} = $split[4] - 1 ; 
		$gene{$id}{"frame"} = $split[7] ;
	}
}
close GFF ; 

### sequence data 
my $chrom = "" ;
my $chr = "start" ;
open FA, "<$fa" ;
while (<FA>) { 
	chomp $_ ; 
	if ( $_ =~ m/^>(.+)/ || eof() ) {
		my $new_chr = $1 ; 
		
		if ( $chr ne "start" ) {
			
            my %syn ;                               ### this hash will store the positions of all 4-d sites on the chromosome
            my %non ;                               ### this hash will store all non-4d sites
			foreach my $mrna ( keys %gene ) {
				if ( $gene{$mrna}{"chrom"} eq $chr ) {
					
                    #### extract the sequences first
					my $protein = "" ;
					my @sites ; 
					if ( $gene{$mrna}{"strand"} eq '+' ) { 
						$protein = substr( $chrom, $gene{$mrna}{"start"}, $gene{$mrna}{"stop"} - $gene{$mrna}{"start"} + 1 ) ; 
						foreach ( $gene{$mrna}{"start"}..$gene{$mrna}{"stop"} ) { 
							push @sites, $_ ;
						}
					}
					else {
                        ### if on - strand, we need to reverse compliment this sequence
						$protein .= reverse_complement ( substr( $chrom, $gene{$mrna}{"start"}, $gene{$mrna}{"stop"} - $gene{$mrna}{"start"} + 1 ) ) ;
						foreach ( $gene{$mrna}{"start"}..$gene{$mrna}{"stop"} ) { 
							push @sites, $_ ;
						}
                        @sites = reverse( @sites ) ; ### reverse list of sites to account for - strand
					}
										
					### identify four fold degenerate sites and non-4d sites
					for ( my $i = $gene{$mrna}{"frame"} ; $i < length( $protein ) - 2; $i += 3 ) { 
						if ( exists( $codon{ substr( $protein, $i, 3 ) } ) ) { 
                            $syn{$chr}{$sites[$i+2]} ++ ;           ### only the third base can be 4d
						}
						else {
							$non{$chr}{$sites[$i+2]} ++ ;
						}
                        $non{$chr}{$sites[$i]} ++ ;                 ### first and second base are always not 4d
						$non{$chr}{$sites[$i+1]} ++ ; 
					}
				}
			}
            
            #### foreach 4d site, check to make sure it's 4d in all transcripts
            #### overlapping reading frames could cause problems here. Hence, whole chromosomes at once
			foreach my $ch ( sort ( keys %syn ) ) { 
				foreach my $pos ( sort { $a <=> $b } ( keys %{$syn{$ch}} ) ) { 
					if ( !exists( $non{$ch}{$pos} ) ) { 
						print $ch, "\t", $pos + 1, "\n" ;
					}
				}
			}
			
            ### reset
			$chrom = "" ;
			$chr = $new_chr ; 
		}
		else { 
			$chr = $new_chr ; 
		}
	}
    
    ### put everything in upper case
	else { 
		$_ =~ s/a/A/g ;
		$_ =~ s/t/T/g ;
		$_ =~ s/c/C/g ;
		$_ =~ s/g/G/g ;
		$chrom .= $_ ; 
	}
}

#### reverse compliment for - strand sequences
sub reverse_complement {
	my $dna = shift;
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}




























