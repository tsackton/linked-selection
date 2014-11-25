#!/bin/bash

for SIZE in 100 500 1000
do
	WINDOW=$(($SIZE*1000))
	./get_pi.pl /Volumes/LaCie/Projects/Current/ne/final_data/poly/4D $WINDOW std
	./get_pi.pl /Volumes/LaCie/Projects/Current/ne/final_data/poly/4D.Q30 $WINDOW q30	
done
