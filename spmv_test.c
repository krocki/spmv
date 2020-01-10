#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "timer.h"
#include "mat.h"
#include "spmat.h"

extern float randf();
extern void randn(float *, float, float, int);
extern void spmat_release(spmat *s);
extern void make_sparse(spmat *spmat, mat *m);
extern void init_mat(mat *m, unsigned M, unsigned K);
extern void spmv(spmat *A, mat *x, mat *b);
extern void mat_print(mat *m, int trans);
extern void sgemm(mat *c, mat *a, mat *b);

void sparsify(mat *m, float rho) {
  for (int i=0; i<((m->ld)*(m->n)); i++) {
    float r = randf();
    if (r > rho) m->d[i] = .0f;
  }
}

unsigned nnz(mat *m) {
  unsigned o=0;
  for (int i=0; i<((m->ld)*(m->n)); i++)
    o += fabs(m->d[i]) > 0.0f;
  return o;
}

#define echo(x...) \
  do { \
    puts(#x); \
    x; \
  } while (0)

struct opts {
  int N, M, K;
  float rho; /* sparsity */
  float max_err;
  int print;
  int bench;
};

int parse_args(int argc, char **argv, struct opts *o);

#define DEFAULT_M 8
#define DEFAULT_N 1
#define DEFAULT_K 8
#define DEFAULT_SPARSITY 0.1f
#define DEFAULT_MAXERR 1e-3f

float cmp(const mat *a, const mat *b) {
  float err;
  float err_norm = .0f;
  for (unsigned i=0; i<a->ld*a->n;i++) {
    err = a->d[i] - b->d[i];
    err *= err;
    err_norm += err;
  }
  return sqrt(err_norm);
}

int main(int argc, char **argv) {

  struct opts opt =  { DEFAULT_N,
                       DEFAULT_M,
                       DEFAULT_K,
                       DEFAULT_SPARSITY,
                       DEFAULT_MAXERR };

  if (0>parse_args(argc, argv, &opt)) return -1;

  int M = opt.M, N = opt.N, K = opt.K;
  float rho = opt.rho;

  double FLOPs = N * M * K * 2.0;
  double mFLOPs = FLOPs / (1 << 20);

  mat a,b,c,f;
  init_mat(&a, M, K);
  init_mat(&b, K, N);
  init_mat(&c, M, N);
  init_mat(&f, M, N);

  randn(a.d, 0, 0.5, M * K);
  randn(b.d, 0, 0.5, K * N);

  double rho_a;
  double rho_b;
  double timer_sgemm;
  double timer_spmv_a;
  double timer_spmv_b;

  if (opt.rho < 1.0f) {
    sparsify(&a, rho);
    sparsify(&b, rho);
  }

  unsigned nnz_a = nnz(&a);
  unsigned nnz_b = nnz(&b);
  rho_a = (double)nnz_a / (double)(M*K);
  rho_b = (double)nnz_b / (double)(K*N);

  int lda = K, ldb = N, ldc = N;

  timeit2(sgemm(&c, &a, &b), &timer_sgemm);

  if (opt.print) {
    echo(mat_print(&a, 1));
    echo(mat_print(&b, 0));
    echo(mat_print(&c, 0));
  }

  spmat sa;
  timeit2(make_sparse(&sa, &a), &timer_spmv_a);
  timeit2(spmv(&sa, &b, &f), &timer_spmv_b);

  float err=cmp(&f, &c);

  if (!opt.bench)
  printf("M %u, K %u, rho %4.3f, T gemm %f, MFLOP/s %.2f, T spmv %f + %f, NNZ A %u, NNZ B %u, err = %f, rho_a %f, rho_b %f\n",
    M, K, rho, timer_sgemm, mFLOPs/timer_sgemm, timer_spmv_b, timer_spmv_a, nnz_a, nnz_b, err, rho_a, rho_b);
  else
  printf("%5u, %5u, %4.3f, %f, %.2f, %f, %f, %5u, %5u, %f, %f, %f\n",
    M, K, rho, timer_sgemm, mFLOPs/timer_sgemm, timer_spmv_b, timer_spmv_a, nnz_a, nnz_b, err, rho_a, rho_b);

  if (opt.print) {
    echo(mat_print(&f, 0));
  }

  free(a.d), free(b.d), free(f.d);
  spmat_release(&sa);
  return 0;
}

int parse_args(int argc, char **argv, struct opts *o) {

  int opt;

  while ((opt = getopt(argc, argv, ":M:N:K:r:hpb")) != -1) {
    switch (opt) {
      case 'M': o->M = strtol(optarg, NULL, 10); break;
      case 'N': o->N = strtol(optarg, NULL, 10); break;
      case 'K': o->K = strtol(optarg, NULL, 10); break;
      case 'r': o->rho = strtof(optarg, NULL); break;
      case 'h': printf("usage: %s -M -N -K -r (sparsity) -b(ench)\n", argv[0]); return -1;
      case 'p': o->print=1; break;
      case 'b': o->bench=1; break;
      default: printf("unknown option %c\n", optopt); return -1;
    }
  }

  for(; optind < argc; optind++) {
    printf("extra arguments: %s\n", argv[optind]);
  }

  return 0;
}
