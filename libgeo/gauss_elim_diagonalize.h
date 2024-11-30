/* gauss_elim_diagonalize.h - Gaussian diagonalization .of a triangular matrix. */
/* Last edited on 2024-11-25 01:19:45 by stolfi */

#ifndef gauss_elim_diagonalize_H
#define gauss_elim_diagonalize_H

#include <stdint.h>

void gauss_elim_diagonalize(uint32_t m, uint32_t n, double M[]);
  /* Assumes that the {m×n} matrix {M} has been totally triangulated,
    Applies row operations to {M} so that it becomes diagonal-like,
    without change to its determinant.
    
    More speciically, let {Lead(M,i)} be the column index of the first
    non-zero element in row {i} of matrix {M}, or {+oo} if that row is
    all zeros. On entry, {Lead(M,i)} must be {+oo} or strictly greater
    than {Lead(M,i-1)} for every {i} (like after a call to
    {gauss_elim_triangularize} with {total=TRUE}). On exit, whenever {j
    := Lead(M,i)} is finite, all elements {M[k,j]} in rows {k < i} are
    set to zero.
    
    The matrix {M} is assumed to be stored as a one-dimensional vector
    {M[0..m*n-1}, in row-by-row order; so that {M[i,j]} is 
    actually {M[i*n + j]}.  */

#endif
