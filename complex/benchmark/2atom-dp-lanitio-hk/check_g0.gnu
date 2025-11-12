
set xrange [0:10]
p '<hgrep latest g0iw-full start 1 1 1 1 1' u 7:8 w lp lt 1, './Greenfunctions/wannier90.green_re' u 1:2 w p lt 3, \
  '<hgrep latest g0iw-full start 1 1 1 1 1' u 7:9 w lp lt 1, './Greenfunctions/wannier90.green_im' u 1:2 w p lt 3, \
  '<hgrep latest g0iw-full start 1 1 1 2 1' u 7:9 w lp lt 1, './Greenfunctions/wannier90.green_im' u 1:3 w p lt 3, \
  '<hgrep latest g0iw-full start 1 1 1 3 1' u 7:9 w lp lt 1, './Greenfunctions/wannier90.green_im' u 1:4 w p lt 3, \
  '<hgrep latest g0iw-full start 1 2 1 5 1' u 7:9 w lp lt 1, './Greenfunctions/wannier90.green_im' u 1:11 w p lt 3, \
  '<hgrep latest g0iw-full start 1 1 1 1 1' u 7:9 w lp lt 1, './Greenfunctions/wannier90.green_im' u 1:237 w p lt 3, \
  '<hgrep latest g0iw-full start 1 2 1 1 1' u 7:8 w lp lt 1, './Greenfunctions/wannier90.green_re' u 1:283 w p lt 3
