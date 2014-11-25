#!/bin/bash

for SIZE in 100 500 1000
do
	WINDOW=$(($SIZE*1000))
	./get_pi_for_windows.pl /Volumes/LaCie/Projects/Current/ne/final_data/vcf/4D $WINDOW std
	./get_pi_for_windows.pl /Volumes/LaCie/Projects/Current/ne/final_data/vcf/4D.Q30 $WINDOW q30	
done
