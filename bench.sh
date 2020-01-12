#!/bin/bash
Ms=(256 4096)
Ks=(256 4096)
As=(0.5 0.05)
Xs=(0.8 0.2)

./spmv_test -l
for k in ${Ks[@]}; do
  for m in ${Ms[@]}; do
    for a in ${As[@]}; do
      for x in ${Xs[@]}; do
        ./spmv_test -M $m -K $k -a $a -x $x -b
      done
    done
  done
done
