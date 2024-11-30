/* gauss_elim_residual.h - Gaussian normalization or a diagonalized matrix. */
/* Last edited on 2024-11-25 01:48:42 by stolfi */

#ifndef gauss_elim_residual_H
#define gauss_elim_residual_H

#include <stdint.h>

void gauss_elim_residual
  ( uint32_t m,
    uint32_t n,
    double A[],
    uint32_t p,
    double B[],
    double X[],
    double R[]
  );
  /* Given the matrices {A} and {B} and a putative solution {X} of the
    system {A X = B}, computes the residual {R = A X - B}. Assumes that
    {R} has size {m × p} (like {B}).
    
    The matrices are assumed to be stored as one-dimensional vectors
    {A[0..m*n-1}, {B[0..m*p-1}, {X[0..n*p-1]}, and {R[0..m*p-1]}, in
    row-by-row order; so that {A[i,j]} is actually {A[i*n + j]}, and
    similarly for the others. */

#endif
