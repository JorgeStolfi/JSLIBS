/* r4x4.h --- 4x4 matrices and operations on them */

#ifndef R4X4_H
#define R4X4_H

#include <stdio.h>
#include "r4.h"

typedef struct { double c[4,4]; } r4x4_t;

void r4x4_map_row (r4_t *a, r4x4_t *m, r4_t *res);
  /* Sets "res" to the product of row vector "a" by matrix "m" */

void r4x4_map_col (r4x4_t *m, r4_t *a, r4_t *res);
  /* Sets "res" to the product of matrix "m" by column vector "a" */

void r4x4_mul (r4x4_t *m, r4x4_t *n, r4x4_t *res);
  /* Sets "res" to the product of matrices "m" and "n" */

double r4x4_cof (r4x4_t *m, int ix, int jx);
  /* Returns the cofactor of element "[ix,jx]" in matrix "m" */

void r4x4_adj (r4x4_t *m, r4x4_t *res);
  /* Sets "res" to the adjoint of matrix "m" */

void r4x4_print (FILE *f, r4x4_t *m);
  /* Prints matrix "m" to file "f" */

#endif

