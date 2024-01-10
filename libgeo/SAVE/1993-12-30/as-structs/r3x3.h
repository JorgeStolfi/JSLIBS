/* r3x3.h --- 3x3 matrices and operations on them */

#ifndef R3X3_H
#define R3X3_H

#include <stdio.h>
#include "r3.h"

typedef struct { double c[3,3]; } r3x3_t;

void r3x3_map_row (r3_t *a, r3x3_t *m, r3_t *res);
  /* Sets "res" to the product of row vector "a" by matrix "m" */

void r3x3_map_col (r3x3_t *m, r3_t *a, r3_t *res);
  /* Sets "res" to the product of matrix "m" by column vector "a" */

void r3x3_mul (r3x3_t *m, r3x3_t *n, r3x3_t *res);
  /* Sets "res" to the product of matrices "m" and "n" */

double r3x3_det (r3x3_t *m);
  /* Returns the determinant of matrix "m" */

double r3x3_cof (r3x3_t *m, int ix, int jx);
  /* Returns the cofactor of element "[ix,jx]" in matrix "m" */

void r3x3_adj (r3x3_t *m, r3x3_t *res);
  /* Sets "res" to the adjoint of matrix "m" */

void r3x3_print (FILE *f, r3x3_t *m);
  /* Prints matrix "m" to file "f" */

#endif

