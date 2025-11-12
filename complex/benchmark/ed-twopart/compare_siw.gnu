#!/usr/local/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 4.6 patchlevel 0    last modified 2012-03-04 
#    	Build System: Linux x86_64
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2012
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal x11  nopersist
# set output
#set term post enh color
#set out 'comare_siw.ps'
set title "Comparison ED vs. CT-Hyb {/Symbol S}(i{/Symbol w})\n  {/Symbol S}(i{/Symbol w}) of the ED calculation is calculated with 1/G_0 - 1/G (the self-energy in sigma.dat is wrong)"
set xlabel "i{/Symbol w}" 
set xrange [ 0 : 5 ] noreverse nowriteback
set ylabel "{/Symbol S}(i{/Symbol w})" 
set yrange [ * : 0 ] noreverse nowriteback
set zero 1e-08
filename=system('ls -ltr *hdf5 | tail -n 1 | sed "s/^.* //"')
p '< hf.py '.filename.' siw  1 1 1 1' u 5:7 w linespo lc rgb '#008000' pt 9 tit 'CT-Hyb {/Symbol S}',\
  '< paste g0mand gm_wim' u 1:(imag($3*{0,1}-1./($5+{0,1}*$6))) w linespo pt 7 lt 3 tit 'ED {/Symbol S}'
#    EOF
