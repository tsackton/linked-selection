Natural selection constrains neutral diversity across a wide range of species
================

This repository contains code to replicate the analysis in "Natural selection constrains neutral diversity across a wide range of species"

Prepare inputs
--------------

This section describes the scripts used to: 
-process and clean genetic maps in preparation for computing genetic maps
-estimate functional density across the genome from GFF files
-estimate theta for four-fold degenerate sites from VCF files

These inputs then form the basis of the population genetic models in the next section.

Initial setup:

1. extract_exons_from_gff.pl: script to convert GFF files with CDS records to bed files containing coding exon records
2. get_gene_density.pl: script to estimate gene density in windows across the genome from exon BED files, launched with run_pi.sh helper script for this analysis
3. get_pi.pl: script to estimate pi in windows across the genome from parsed and filtered VCF files, launched with run_gd.sh helper script for this analysis
4. load_maps.R: script to clean up genetic maps and produce input for recombination rate estimation; run as R --vanillia < load_maps.R

Recombination rate estimation
-----------------------------





Population genetics models
--------------------------

This section describes the scripts used to:
-compute the effect of background selection on windows across the genome, taking into account recombination rate and gene density,
based on the method of Rockman et al and Flowers et al.
-fit a variety of population genetic models to the data using nonlinear regression in R


Analysis
--------

This section describes the scripts and code used to:
-Process the output of the population genetics models above
-Estimate species range from occurrence data
-Replicate the analyses presented in the manuscript
