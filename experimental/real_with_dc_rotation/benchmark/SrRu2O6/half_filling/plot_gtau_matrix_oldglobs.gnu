#!/usr/bin/gnuplot --persist
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
set yrange [0.0001:1]
set xrange [  ] noreverse nowriteback
#set logscale y
p "Gt_QMC_first_iteration.dat" u 1:(-$2) w l lc rgb "black",\
  "first_iteration/gtau_matrix_oldglobs.dat" u 5:6 index 0 w p  lc rgb 'red' ps 1,\
  "first_iteration/gtau_matrix_oldglobs.dat" u 5:6 index 2 w p  lc rgb 'blue' ps 1,\
  "first_iteration/gtau_matrix_oldglobs.dat" u 5:6 index 4 w p  lc rgb 'dark-olivegreen' ps 1
