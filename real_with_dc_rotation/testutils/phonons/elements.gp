set term pngcairo size 800,800
set out "colormatrix.png"
set yrange [:] reverse
set palette rgbformula 2,-7,-7
set logscale cb
set size square
unset colorbox
set autoscale fix
set tics scale 0
unset cbtics
set cblabel 'Score'
set xtics format ""

set arrow 1 from 0.5,4.5 to 12.5,4.5 nohead front
set arrow 2 from 4.5,0.5 to 4.5,12.5 nohead front

set arrow 3 from 2.5,0.5 to 2.5,4.5 nohead front dt 2
set arrow 4 from 0.5,2.5 to 4.5,2.5 nohead front dt 2


set arrow 5 from 8.5,0.5 to 8.5,12.5 nohead front dt 2
set arrow 6 from 0.5,8.5 to 12.5,8.5 nohead front dt 2

set x2tics (1,2,3,4,5,6,7,8,9,10,11,12)
set ytics (1,2,3,4,5,6,7,8,9,10,11,12)
unset key
plot 'elements.dat' u ($1+1):($2+1):(abs($3)) matrix with image,\
     '' matrix using ($1+1):($2+1):(sprintf('%.2f', $3)) with labels font ',11'
