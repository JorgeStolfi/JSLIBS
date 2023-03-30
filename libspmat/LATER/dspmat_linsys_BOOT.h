#ifndef dspmat_linsys_BOOT_H
#define dspmat_linsys_BOOT_H
/* Gauss-Seidel linear system solving for sparse matrices with {double} entries */

#define dspmat_linsys_BOOT_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-08-17 by J.Stolfi, UNICAMP */
/* Last edited on 2023-03-18 10:51:12 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <dspmat.h>

/* !!! BUGGY -- need preconditioning of the {M} matrix. !!! */

void dspmat_linsys_BOOT_solve
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
  /* Computes the solution {x} of the linear equation system {b = A
    x}, where {A} is a given {nb × nx} sparse matrix, {b[0..nb-1]} is
    a given (column) vector, and {x[0..nx-1]} is an unknown (column)
    vector.
    
    Requires the matrix to be square, hence {nb == nx}. Said another
    way, stores into {x} the product {A^{-1}*b}.  The matrix must not
    contain two entries with the same indices.
    
    Uses a bootstrap method method with parameters
    {omega,abs_tol,rel_tol}, starting with {x = (0,0,.. 0)}. Stops
    after doing {max_iter} passes, or as soon as the residual is small
    according to {abs_tol,rel_tol}. */ 

void dspmat_inv_mul_BOOT
  ( dspmat_t *A,
    dspmat_t *B,
    dspmat_t *X,
    int32_t max_iter, 
    double abs_tol, 
    double rel_tol
  );
  /* Sets {*X} to {A^{-1}*B}, the inverse of {A} times the matrix {B}; namely,
    the solution of the matrix equation {A*X == B}. The matrix {A} must be
    square. 
    
    Uses an interative method; stops after {max_iter} iterations, or
    when the elements of the residual {B - A*X} are small according to
    {abs_tol,rel_tol}. */

#endif
