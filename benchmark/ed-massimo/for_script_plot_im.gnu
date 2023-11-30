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
#p "fort.31" u 1:($1*$3) w p lc rgb "grey" pt 7,\
 #"fort.32" u 1:($1*$3) w p lc rgb "black" pt 7,\
#p "siw111_beta100.dat" u 1:3 w p lc rgb "grey" pt 7,\
 #"siw121_beta100.dat" u 1:3 w p lc rgb "black" pt 7,\
 #'< hgrep latest siw 1 1 1 1' u 5:($5*$7) w l lc rgb "red",\
 #'< hgrep latest siw 1 1 2 1' u 5:($5*$7) w l lc rgb "blue"
set term png enhanced size 1280,1024
set output "out.png"
set yrange [-3:0]
set xrange [0:50]
p -1.47 lc rgb "grey" pt 7,\
 -1.39 lc rgb "black" pt 7,\
 '< hgrep latest siw 1 1 1 1' u 5:($5*$7) w l lc rgb "red",\
 '< hgrep latest siw 1 1 2 1' u 5:($5*$7) w l lc rgb "blue"
