#GNUPLOT
set key bott
set xrange [0:30]
p'< hgrep benchmark.hdf5 siw -1 1 1 1' u 5:6 w lp, 'CTINT/Sigma_iom_qmc' u 1:($2+3.34/2.) w lp lt 3
