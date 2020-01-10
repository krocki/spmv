#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mat.h"

void spmv_nn(
 const int M, const int K,
 float *restrict a, const int lda,
 float *restrict b, const int ldb,
 float *restrict c, const int ldc) {

  for (int i=0; i<M; i++)
  for (int k=0; k<K; k++)
    c[i*ldc] += a[i*lda+k] * b[k*ldb];

}

void sgemm_nn(
 const int M, const int N, const int P,
 float *restrict a, const int lda,
 float *restrict b, const int ldb,
 float *restrict c, const int ldc) {

  for (int i=0; i<M; i++)
  for (int k=0; k<P; k++)
  for (int j=0; j<N; j++)
    c[i*ldc+j] +=
    a[i*lda+k] * b[k*ldb+j];

}

void sgemm(mat *c, mat *a, mat *b) {

  int at = 0;
  int bt = 0;
  int M = c->n;
  int N = c->ld;
  int K = a->ld;
  if (1==N)
    spmv_nn(M, K, a->d, a->ld, b->d, b->ld, c->d, c->ld);
  else
    sgemm_nn(M, N, K, a->d, a->ld, b->d, b->ld, c->d, c->ld);
  //cblas
  //int at = a->transp;
  //int bt = b->transp;

  //unsigned M = c->rows;
  //unsigned N = c->cols;
  //unsigned K = at ? a->rows : a->cols;

  //float alpha = 1.0f;
  //float beta = 1.0f;

  //unsigned lda = (1==at) ? K : M;
  //unsigned ldb = (1==bt) ? N : K;
  //unsigned ldc = M;

  //cblas_sgemm( col_major, at, bt, M, N, K, alpha,
  //             a->d, lda, b->data, ldb,
  //             beta, c->data, ldc);

}
