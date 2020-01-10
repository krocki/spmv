# SpMV
Sparse Matrix-Vector (SpMV) bench.

`x` is a dense array

`y` is a dense array

`A` is a sparse matrix.

```y = Ax```

```
./spmv -M 1024 -K 1024 -r 0.05
```

```
M 1024, K 1024, rho 0.050, T gemm 0.026158, MFLOP/s 76.46, T spmv 0.000018 + 0.013132, NNZ A 52465, NNZ B 44, err = 0.000000, rho_a 0.050035, rho_b 0.042969
```
