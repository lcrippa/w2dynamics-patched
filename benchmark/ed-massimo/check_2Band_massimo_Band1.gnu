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
set xrange [ 0.00000 : 5.00000 ] noreverse nowriteback
filename=system('ls -ltr *hdf5 | tail -n 1 | sed "s/^.* //"')
p '< hgrep '.filename.' siw 1 1 1 1' u 5:7 w linespo lc rgb '#008000' pt 9,\
  'siw111_beta100.dat' u 1:3 w linespo pt 7 lt 3
#    EOF
