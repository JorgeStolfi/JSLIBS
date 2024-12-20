/* rmxn_spin.h --- applying random rotations to matrices */
/* Last edited on 2024-12-05 10:28:39 by stolfi */

#ifndef rmxn_spin_H
#define rmxn_spin_H

#include <stdint.h>

void rmxn_spin_rows(uint32_t m, uint32_t n, double A[], double M[]);
  /* Applies a random rotation to each row of {A}, which is assumed to
    have {m} rows and {n} columns. Equivalent to computing {M=A*N}
    where {N} is a random orthonormal {n×n} matrix. */

void rmxn_spin_cols(uint32_t m, uint32_t n, double A[], double M[]);
  /* Applies a random rotation to each column {A}, which is assumed to
    have {m} rows and {n} columns. Equivalent to computing {M=N*A}
    where {N} is a random orthonormal {m×m} matrix. */

#endif

