#!/usr/bin/gnuplot -persist
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
p \
"../norot/siw_rot00.dat" u 1:3 w l lc rgb "black",\
"../norot/siw_rot01.dat" u 1:3 w l lc rgb "black",\
"../norot/siw_rot10.dat" u 1:3 w l lc rgb "black",\
"../norot/siw_rot11.dat" u 1:3 w l lc rgb "black",\
"../norot/siw_rot22.dat" u 1:3 w l lc rgb "black",\
"../norot/siw_rot23.dat" u 1:3 w l lc rgb "black",\
"../norot/siw_rot32.dat" u 1:3 w l lc rgb "black",\
"../norot/siw_rot33.dat" u 1:3 w l lc rgb "black",\
"siw-full.dat" u 7:($9) w p pt 7 ps 0.5 lc rgb "blue"
