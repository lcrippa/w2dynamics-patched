#!/usr/bin/gnuplot -persist
#
#    
#     G N U P L O T
#     Version 4.6 patchlevel 0    last modified 2012-03-04 
#     Build System: Linux x86_64
#    
#     Copyright (C) 1986-1993, 1998, 2004, 2007-2012
#     Thomas Williams, Colin Kelley and many others
#    
#     gnuplot home:     http://www.gnuplot.info
#     faq, bugs, etc:   type "help FAQ"
#     immediate help:   type "help"  (plot window: hit 'h')
# set terminal x11  nopersist
set yrange [-1:1]
set xrange [-100:100]
p \
'< hgrep sc-2017-02-23-Thu-14-54-46.hdf5 siw-full -2 1 1 1 1 1 ' u 7:9 w l lc rgb "grey",\
'< hgrep sc-2017-02-23-Thu-14-54-46.hdf5 siw-full -2 1 1 2 1 2 ' u 7:9 w l lc rgb "grey",\
'< hgrep sc-2017-02-23-Thu-14-54-46.hdf5 siw-full -2 1 1 1 1 2 ' u 7:9 w l lc rgb "grey",\
'< hgrep sc-2017-02-23-Thu-14-54-46.hdf5 siw-full -2 1 1 2 1 1 ' u 7:9 w l lc rgb "grey",\
'< hgrep sc-2017-02-23-Thu-14-54-46.hdf5 siw-full -1 1 1 1 1 1 ' u 7:9 w l lc rgb "grey",\
'< hgrep sc-2017-02-23-Thu-14-54-46.hdf5 siw-full -1 1 1 2 1 2 ' u 7:9 w l lc rgb "grey",\
'< hgrep sc-2017-02-23-Thu-14-54-46.hdf5 siw-full -1 1 1 1 1 2 ' u 7:9 w l lc rgb "grey",\
'< hgrep sc-2017-02-23-Thu-14-54-46.hdf5 siw-full -1 1 1 2 1 1 ' u 7:9 w l lc rgb "grey",\
"../norot/siw_rot00.dat" u 1:3 w p pt 7 ps 0.5 lc rgb "blue",\
"../norot/siw_rot01.dat" u 1:3 w p pt 7 ps 0.5 lc rgb "blue",\
"../norot/siw_rot10.dat" u 1:3 w p pt 7 ps 0.5 lc rgb "blue",\
"../norot/siw_rot11.dat" u 1:3 w p pt 7 ps 0.5 lc rgb "blue"
