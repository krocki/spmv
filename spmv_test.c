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
extern unsigned long spmv_flops(spmat *A, mat *x);
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
  float rho_a; /* A sparsity */
  float rho_x; /* x sparsity */
  float max_err;
  int print;
  int bench;
  int legend;
};

int parse_args(int argc, char **argv, struct opts *o);

#define DEFAULT_M 256
#define DEFAULT_N 1
#define DEFAULT_K 512
#define DEFAULT_SPARSITY 0.5f
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
                       DEFAULT_SPARSITY,
                       DEFAULT_MAXERR };

  if (0>parse_args(argc, argv, &opt)) return -1;

  int M = opt.M, N = opt.N, K = opt.K;
  float rho_a_target = opt.rho_a;
  float rho_x_target = opt.rho_x;

  mat a,b,c,f;
  init_mat(&a, M, K);
  init_mat(&b, K, N);
  init_mat(&c, M, N);
  init_mat(&f, M, N);

  randn(a.d, 0, 0.5, M * K);
  randn(b.d, 0, 0.5, K * N);

  double rho_a;
  double rho_x;
  double timer_sgemm;
  double timer_spmv_a;
  double timer_spmv_b;

  if (rho_a_target < 1.0f)
    sparsify(&a, rho_a_target);
  if (rho_x_target < 1.0f)
    sparsify(&b, rho_x_target);

  unsigned nnz_a = nnz(&a);
  unsigned nnz_b = nnz(&b);
  double dFLOPs = N * M * K * 2.0;
  double mFLOPs = dFLOPs / (1 << 20);

  rho_a = (double)nnz_a / (double)(M*K);
  rho_x = (double)nnz_b / (double)(K*N);

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

  double sp_flops = spmv_flops(&sa, &b);
  sp_flops /= (1 << 20);
  float err=cmp(&f, &c);

  if (!opt.bench)
  printf("M %u, K %u, ra %4.3f, rx %4.3f\ngemm/spmv time %2.5f, T gemm %f\ngemm_mflops = %4.6f, d MFLOP/s %.2f, T spmv %f + %f\nspmv_mflops = %4.6f, s MFLOP/s %.2f, NNZ A %u, NNZ B %u\nrho_a %f, rho_x %f, err = %f\n", M, K, rho_a_target, rho_x_target, timer_sgemm / timer_spmv_b, timer_sgemm, mFLOPs, mFLOPs/timer_sgemm, timer_spmv_b, timer_spmv_a, sp_flops, sp_flops / timer_spmv_b, nnz_a, nnz_b, rho_a, rho_x, err);
  else
  printf("%5u, %5u, %3.2f, %3.2f, %7.3f, %7f, %4.1f, %7.1f, %8f, %7f, %5.2f, %7.1f, %8u, %6u, %4.2f, %4.2f, %f\n", M, K, rho_a_target, rho_x_target, timer_sgemm/timer_spmv_b, timer_sgemm, mFLOPs, mFLOPs/timer_sgemm, timer_spmv_b, timer_spmv_a, sp_flops, sp_flops/timer_spmv_b, nnz_a, nnz_b, rho_a, rho_x, err);

  if (opt.print) {
    echo(mat_print(&f, 0));
  }

  free(a.d), free(b.d), free(f.d);
  spmat_release(&sa);
  return 0;
}

void print_tab_legend() {
  printf("%5s, %5s, %4s, %4s, %7s, %8s, %4s, %4s, %3s, %7s, %5s, %7s, %8s, %6s, %4s, %4s, %s\n",
  "M" , "K" , "ra" , "rx" , "X" , "t sgemm" , "mF" , "gemm FL" , "t spmv b" , "t spmv a" , "ops" , "sp mFL" , "nnz_a" , "nnz_b" , "r_a" , "r_x" , "err");
}

int parse_args(int argc, char **argv, struct opts *o) {

  int opt;

  while ((opt = getopt(argc, argv, ":M:N:K:a:x:hpbl")) != -1) {
    switch (opt) {
      case 'M': o->M = strtol(optarg, NULL, 10); break;
      case 'N': o->N = strtol(optarg, NULL, 10); break;
      case 'K': o->K = strtol(optarg, NULL, 10); break;
      case 'a': o->rho_a = strtof(optarg, NULL); break;
      case 'x': o->rho_x = strtof(optarg, NULL); break;
      case 'h': printf("usage: %s -M -N -K -a (sparsity) -x (sparsity) -b(ench)\n", argv[0]); return -1;
      case 'p': o->print=1; break;
      case 'l': print_tab_legend(); return -1;
      case 'b': o->bench=1; break;
      default: printf("unknown option %c\n", optopt); return -1;
    }
  }

  for(; optind < argc; optind++) {
    printf("extra arguments: %s\n", argv[optind]);
  }

  return 0;
}
