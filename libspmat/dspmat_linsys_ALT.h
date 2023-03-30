#ifndef dspmat_linsys_ALT_H
#define dspmat_linsys_ALT_H
/* Iterative linear system solving for sparse matrices with {double} entries */

#define dspmat_linsys_ALT_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-08-17 by J.Stolfi, UNICAMP */
/* Last edited on 2023-03-18 10:49:34 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>
#include <dspmat.h>

void dspmat_linsys_ALT_solve
  ( double b[],
    dspmat_size_t nb,
    dspmat_t *A, 
    double x[], 
    dspmat_size_t nx,
    int32_t max_iter, 
    double omega,
    double abs_tol, 
    double rel_tol
  );
  /* Same as {dspmat_linsys_solve} but using an alternative algorithm. */
    
double dspmat_linsys_ALT_compute_residuals
  ( double b[],
    dspmat_t *A, 
    double x[],
    double d[],
    double abs_tol, 
    double rel_tol
  );
  /* Sets {d[i]} to the residual {b[i] - SUM{A[i,j]*x[j] : j}}, for
    each {i} in {0..n-1} where {n == A->rows == A->cols}. Returns the
    magnitude {err} of the residual relative to the absolute and
    relative tolerances {abs_tol,rel_tol} (so that {err <= 1.0} means
    `within tolerance'). */

void dspmat_linsys_ALT_recompute_solution
  ( double b[],
    dspmat_t *A,
    double x[], 
    double d[],
    double c[],
    double r[],
    double M[]
  );
  /* Recomputes the solution {x[0..n-1]} from the the current guess
    {x[0..n-1]} and the coresponding residuals {d[0..n-1]}, where {n
    == A->rows == A->cols}. Requires the col and row sums {c[0..n-1],
    r[0..n-1]} of {A}, and a work area {M[0..ne-1]} where {ne ==
    A->ents}. */

#endif
