/* r2x2.h --- 2x2 matrices and operations on them */
/* Last edited on 2001-10-21 21:17:48 by stolfi */

#ifndef r2x2_H
#define r2x2_H

#include <stdio.h>
#include "r2.h"

typedef double r2x2_t [2][2];

void r2x2_map_row (r2_t a, r2x2_t m, r2_t res);
  /* Sets "res" to the product of row vector "a" by matrix "m" */

void r2x2_map_col (r2x2_t m, r2_t a, r2_t res);
  /* Sets "res" to the product of matrix "m" by column vector "a" */

void r2x2_mul (r2x2_t m, r2x2_t n, r2x2_t res);
  /* Sets "res" to the product of matrices "m" and "n" */

double r2x2_det (r2x2_t m);
  /* Returns the determinant of matrix "m" */

double r2x2_cof (r2x2_t m, int ix, int jx);
  /* Returns the cofactor of element "[ix,jx]" in matrix "m" */

void r2x2_adj (r2x2_t m, r2x2_t res);
  /* Sets "res" to the adjoint of matrix "m" */

void r2x2_print (FILE *f, r2x2_t m);
  /* Prints matrix "m" to file "f" */

#endif

