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
set yrange [-1:1]
set xrange [-100:100]
p "atsushi/Sigma.dat" u 1:3 index 0 w l lc rgb "black",\
  "atsushi/Sigma.dat" u 1:3 index 1 w l lc rgb "grey",\
 "siw_segment.dat" index 0 u 5:7 w p ps 1 ,\
 "siw_segment.dat" index 1 u 5:7 w p ps 1 ,\
 "siw_segment.dat" index 2 u 5:7 w p ps 1 
#"siw_matrix.dat" index 0 u 5:7 w p ps 1 ,\
#"siw_matrix.dat" index 1 u 5:7 w p ps 1 ,\
#"siw_matrix.dat" index 2 u 5:7 w p ps 1 
