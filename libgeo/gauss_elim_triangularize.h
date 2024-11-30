/* gauss_elim_triangularize.h - Gaussian triangulation of a matrix. */
/* Last edited on 2024-11-25 01:42:24 by stolfi */

#ifndef gauss_elim_triangularize_H
#define gauss_elim_triangularize_H

#include <stdint.h>

#include <bool.h>

void gauss_elim_triangularize
  ( uint32_t m,
    uint32_t n,
    double M[],
    bool_t total,
    double tiny
  );
  /* Applies the Gaussian elimination method to the {m×n} matrix
    {M}, leaving it upper triangular. 
    
    Specifically, the matrix {M} is modified by row operations that do
    not change its determinant. Let {Lead(M,i)} be the column index of
    the first non-zero element in row {i} of matrix {M}, or {+oo} if
    that row is all zeros. Upon exit, in any case, we will have
    {Lead(M,i) >= i} for every {i}. Thus, in particular, if {m = n},
    the determinant of {M} will be the product of its diagonal
    elements {M[i,i]} for {i} in {0..m-1}.
    
    If {total} is true, the output matrix satisfies a stronger
    condition: for all {i} in {1..m-1}, {Lead(M,i) < Lead(M,i-1)}, or
    {Lead(M,i) == Lead(M,k) = +oo}. In this case, the non-zero rows of
    {M} are a basis for the space spanned by the original rows.
    
    During triangularization, any entry whose absolute value gets
    reduced to {tiny} or less times its previous value is set to zero.
    If {tiny} is zero or negative, this cleanup is supressed.
    
    The matrix {M} is assumed to be stored as a one-dimensional vector
    {M[0..m*n-1}, in row-by-row order; so that {M[i,j]} is 
    actually {M[i*n + j]}. */

#endif
