#GNUPLOT
set key top
set xrange [ 0.00000 : 20.0000 ] noreverse nowriteback
set yrange [*:0]
p'< hgrep benchmark.hdf5 siw -1 1 1 1' u 5:7 w lp, 'CTINT/Sigma_iom_qmc' u 1:3 w lp lt 3
