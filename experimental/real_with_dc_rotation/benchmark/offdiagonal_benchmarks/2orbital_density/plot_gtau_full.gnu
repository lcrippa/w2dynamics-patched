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
#asdf
#set terminal postscript eps color enhanced  
#set output "G_full.eps"
#set ylabel "G({/Symbol t})"
#set xlabel "{/Symbol t}"
#set key bottom center
p \
"ED/ed-result/gtau.dat" u 3:(-$4) w p pt 7 ps 0.5 lc rgb "black",\
'< hgrep latest gtau-full 1 1 1 1 1 1 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 1 1 1 2 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 1 1 2 1 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 1 1 2 2 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 1 2 1 1 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 1 2 1 2 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 1 2 2 1 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 1 2 2 2 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 2 1 1 1 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 2 1 1 2 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 2 1 2 1 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 2 1 2 2 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 2 2 1 1 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 2 2 1 2 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 2 2 2 1 ' u 7:8 w l ,\
'< hgrep latest gtau-full 1 1 2 2 2 2 ' u 7:8 w l 
