#!/bin/bash
Ms=(4 16 64 128 1024)
Ks=(4 16 64 128 1024 2048)
Rs=(0.01 0.05 0.1 0.2 0.5 1.0)

for m in ${Ms[@]}; do
  for k in ${Ks[@]}; do
    for r in ${Rs[@]}; do
      ./spmv_test -M $m -K $k -r $r -b
    done
  done
done
