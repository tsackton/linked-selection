#!/bin/bash

for SIZE in 100 500 1000
do
	WINDOW=$(($SIZE*1000))
	./get_gene_density_for_windows.pl /Volumes/LaCie/Projects/Current/ne/final_data/gd $WINDOW
done