#!/bin/bash
Ms=(32 64 96 128 256 512 768 1024 2048 4096 8192)
Ks=(32 64 96 128 256 512 758 1024 2048 4096 8192)
As=(1.0 0.5 0.2 0.1 0.01 0.02 0.01 0.005 0.001)
Xs=(1.0 0.75 0.5 0.25 0.1 0.05 0.01)

./spmv_test -l
for a in ${As[@]}; do
  for x in ${Xs[@]}; do
    for m in ${Ms[@]}; do
      for k in ${Ks[@]}; do
        ./spmv_test -M $m -K $k -a $a -x $x -b
      done
    done
  done
done
