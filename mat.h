#ifndef MAT_H
#define MAT_H

typedef struct {
  float *d;
  int ld;
  int n;
} mat;

#define IJ(m,y,x) &((m).d[((y)*((m).ld)+(x))])

#endif
