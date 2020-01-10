#!/bin/bash
Ms=(64 128 1024 2048 4096 8192)
Ks=(64 128 1024 2048 4096 8192)
As=(1.0 0.2 0.05 0.01)
Xs=(1.0 0.5 0.2 0.1)

for a in ${As[@]}; do
  for x in ${Xs[@]}; do
    for m in ${Ms[@]}; do
      for k in ${Ks[@]}; do
        ./spmv_test -M $m -K $k -a $a -x $x -b
      done
    done
  done
done
