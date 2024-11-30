/* gauss_elim_normalize.h - Gaussian normalization or a diagonalized matrix. */
/* Last edited on 2024-11-25 01:20:27 by stolfi */

#ifndef gauss_elim_normalize_H
#define gauss_elim_normalize_H

#include <stdint.h>

void gauss_elim_normalize(uint32_t m, uint32_t n, double M[]);
  /* Assumes that the {m×n} matrix {M} has been totally triangularized. 
    Scales each row so that the leading nonzero element is 1.0.
    
    More speciically, let {Lead(M,i)} be the column index of the first
    non-zero element in row {i} of matrix {M}, or {+oo} if that row is
    all zeros. On entry, {Lead(M,i)} must be {+oo} or strictly greater
    than {Lead(M,i-1)} for every {i} (like after a call to
    {gauss_elim_triangularize} with {total=TRUE}).
    
    On exit, whenever {j:=Lead(M,i)} is finite, scales row {i} so that {M[i,j] == 1.0}.
    
    The matrix {M} is assumed to be stored as a one-dimensional vector
    {M[0..m*n-1}, in row-by-row order; so that {M[i,j]} is 
    actually {M[i*n + j]}. */

#endif
