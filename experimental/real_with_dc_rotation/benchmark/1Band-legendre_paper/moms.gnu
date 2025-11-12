#GNUPLOT
set key bott
set xrange [ 0.00000 : 50.0000 ] 
set logsc y
set yrange [ 0.000100000 : 100.000 ] 
D = 2.0
u = 4.0
mom = 5.0
p 'c.dat' u ($0*2):1 pt 7 ps 2, 1 lt 0 , 'cc.dat' u ($0*2+1):1 pt 7 ps 2 lt 2, 'ccc.dat' u ($0*2):1 pt 7 ps 2 lt 3, mom lt 0 lw 4, 0 lt 0
