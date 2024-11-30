/* gauss_elim_determinant.h - Gaussian normalization or a diagonalized matrix. */
/* Last edited on 2024-11-25 01:46:02 by stolfi */

#ifndef gauss_elim_determinant_H
#define gauss_elim_determinant_H

#include <stdint.h>

double gauss_elim_determinant(uint32_t m, uint32_t n, double A[], uint32_t q);
  /* Returns the determinant of the first {q} rows and {q} columns of
    the {m × n} array {A}. If {q > m} or {q > n}, the result is zero.
    The array {A} is not changed. */

double gauss_elim_determinant_triang(uint32_t m, uint32_t n, double M[], uint32_t q);
  /* Assumes that the {m×n} matrix {M} has been triangularized
    as per {gauss_elim_triangularize}, with {total} either true or false.
    Returns the product of the elements on the main diagonal of the
    first {q} rows and columns. Returns zero if {q > m} or {q > n}.
    Note that when {q == m} the result is the determinant of the first
    {m} columns of {M}.
    
    The matrix {M} is assumed to be stored as a one-dimensional vector
    {M[0..m*n-1}, in row-by-row order; so that {M[i,j]} is 
    actually {M[i*n + j]}. */

#endif
