#!/bin/sh

if [ -z $maxg ] ; then
	echo "plz input maxg"
	exit 1
fi

if [ -z $filename ] ; then
	echo "plz input filename"
	exit 1
fi

today=$(date "+%y%m%d_%H%M")

gnuplot -e "
	unset key ;
	set terminal png ;
	set output 'Plot_Obj_${today}.png' ;
	plot [0:${maxg}] '${filename}' using 1:2 with lines lc 'black'
	"

gnuplot -e "
	unset key ;
	set terminal png ;
	set output 'Plot_Peri_${today}.png' ;
	plot [0:${maxg}] '${filename}' using 1:3 with lines lc 'black'
	"

gnuplot -e "
	unset key ;
	set terminal png ;
	set output 'Plot_Fit_${today}.png' ;
	plot [0:${maxg}] '${filename}' using 1:3 with lines lc 'black'
	"

