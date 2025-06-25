#GNUPLOT
#set term post enh color solid "Times-Roman" 24
#set out "sigma_U5.ps"
set title "U=5 {/Symbol b}=50"
set xlabel '{/Symbol w}_n'
set ylabel 'Im {/Symbol S}({/Symbol w}_n)'
set key bott
set xrange [ 0 : 2 ]
set yrange [ * : 0.00000 ] 
p 'beta_50/self-en_wim' u 1:3 w l lw 3 lt 3 title "ED", \
  '< hgrep latest siw -1 1 1 1' u 5:7 w p pt 7 ps 1.5 lt 1  title "CT"
