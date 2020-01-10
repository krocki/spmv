#ifndef SP_MAT_H
#define SP_MAT_H

typedef struct {
  unsigned i;
  float w;
} edge;

typedef struct {
  edge *e;
  unsigned len;
  unsigned cap;
} vec;

typedef struct {
  vec *cols;
  unsigned n_cols;
} spmat;

#endif
