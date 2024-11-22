/* rmxn_shift.h --- adding scalars to rows or cols of MxN matrices */
/* Last edited on 2024-11-20 14:08:33 by stolfi */

#ifndef rmxn_shift_H
#define rmxn_shift_H

#define _GNU_SOURCE
#include <stdint.h>

/* OTHER */

void rmxn_shift_rows(uint32_t m, uint32_t n, double A[], double v[], double M[]);
  /* Adds the vector {v[0..n-1]} to each row of the matrix {A},
    yielding matrix {M}. Both matrices are assumed to have {m} rows
    and {n} columns. */

void rmxn_shift_cols(uint32_t m, uint32_t n, double A[], double v[], double M[]);
  /* Adds the vector {v[0..m-1]} to each column of the matrix {A},
    yielding matrix {M}. Both matrices are assumed to have {m} rows
    and {n} columns. */

#endif

