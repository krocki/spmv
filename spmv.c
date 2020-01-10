#include "spmat.h"
#include "mat.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/* helper vector functions */
void vec_init(vec *v, unsigned cap) {
  v->e = calloc(cap, sizeof(edge));
  if (NULL != v->e) {
    v->cap=cap;
    v->len=0;
  }
}

void vec_release(vec *v) {
  if (NULL != v->e) {
    free(v->e); v->e = NULL;
    v->cap = 0; v->len = 0;
  }
}

void vec_resize(vec *v, unsigned cap) {
  v->e = realloc(v->e, sizeof(edge) * cap);
  v->cap = cap;
}

void vec_append(vec *v, edge *e) {
  if ((v->len+1) >= v->cap) {
    vec_resize(v,
      v->cap == 0 ? 1 :
      2 * v->cap);
  }
  v->e[v->len++] = *e;
}

/* b := Ax */
void spmv(spmat *A, mat *x, mat *b) {

  for (unsigned i=0; i<x->n; i++) {
    float x_val = x->d[i];
    if (fabs(x_val) > 0) {
      vec *a_col = &(A->cols[i]);
      unsigned nnz_a = a_col->len;
      for (unsigned j=0; j<nnz_a; j++) {
        edge *e = &(a_col->e[j]);
        float a_val = e->w;
        unsigned a_i = e->i;
        b->d[a_i] += x_val * a_val;
      }
    }
  }
}

unsigned long spmv_flops(spmat *A, mat *x) {
  unsigned long flops = 0UL;
  for (unsigned i=0; i<x->n; i++) {
    float x_val = x->d[i];
    if (fabs(x_val) > 0) {
      vec *a_col = &(A->cols[i]);
      unsigned nnz_a = a_col->len;
      for (unsigned j=0; j<nnz_a; j++) {
        flops += 2;
      }
    }
  }
  return flops;
}

void print_sparse(spmat *s) {
  for (unsigned i=0; i<s->n_cols; i++) {
    if (s->cols[i].len>0) printf("col %u: ", i);
    for (unsigned j=0; j<s->cols[i].len; j++) {
      edge *e = &(s->cols[i].e[j]);
      printf("(%u, %.3f)%s",
        e->i, e->w, (j==(s->cols[i].len-1)) ? "\n":" ");
    }
  }
}

void init_mat(mat *m, unsigned M, unsigned K) {

  m->d = calloc(M * K, sizeof(float));
  m->ld = K;
  m->n = M;

}

void make_sparse(spmat *spmat, mat *m) {

  spmat->n_cols = m->ld;
  spmat->cols = (vec*)calloc(m->ld,
    sizeof(vec));

  for (unsigned c=0; c<m->ld; c++) {
    for (unsigned r=0; r<m->n; r++) {
      float v = m->d[r * m->ld + c];
      if (fabs(v) > .0f) {
        edge e = {r, v};
        vec_append(&(spmat->cols[c]), &e);
      }
    }
  }
}

void spmat_release(spmat *s) {
  for (unsigned i=0; i<s->n_cols; i++)
    vec_release(&(s->cols[i]));
  free(s->cols);
}

void arr_print(float *arr, size_t rows, size_t cols, size_t ld) {
  for (size_t i=0; i<rows; i++) {
    for (size_t j=0; j<cols; j++) {
      printf("%8.5f%c",
       arr[j+i*ld], j==(cols-1) ? '\n' : ' ');
    }
  }
}

void mat_print(mat *m, int trans) {
  arr_print(m->d, trans ? m->n : m->ld, trans ? m->ld : m->n, m->ld);
}
