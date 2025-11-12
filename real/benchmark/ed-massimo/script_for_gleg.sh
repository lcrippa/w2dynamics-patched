


nleg_min=1
nleg_max=100


for((i=$nleg_min; i<$nleg_max+1; i++))
do
   echo " ========================================= durchlauf $i ==============================================="
   DMFT.py Parameters.in $i
   gnuplot --persist for_script_plot_im.gnu
   #filename=`ls -ltr | tail -n 1 | sed "s/^.* //"`
   #echo $filename
   k="$(printf '%03d' "$i")"
   #mv out.png out_$k.png
   convert out.png out_$k.pdf
   rm out.png
done

#for name in out*; do convert $name `basename $name .png`.pdf; done
#rm *.png

convert out* output.pdf
