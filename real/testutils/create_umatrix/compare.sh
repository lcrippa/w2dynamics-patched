#!/bin/zsh

NLINES=$(cat u_matrix_template.dat | wc -l)

for i in {3..$NLINES}; do
  LINE=$(head -n $i u_matrix_template.dat | tail -n 1)
  INDICES=$(echo $LINE | awk '{print $1" "$2" "$3" "$4" "$5" "$6" "$7" "$8}')
  MATCHING=$(cat u_matrix.dat | grep $INDICES | wc -l)
  if [ $MATCHING = "0" ]; then
    echo "LINE $i DOESN'T MATCH"
    exit
  fi
done
echo "FILES MATCH!"

