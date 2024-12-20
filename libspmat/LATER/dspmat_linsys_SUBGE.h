#ifndef dspmat_linsys_SUBGE_H
#define dspmat_linsys_SUBGE_H
/* Linear system solving by Sub-Gaussian-elimination for sparse {double} matrices. */

#define dspmat_linsys_SUBGE_H_COPYRIGHT "Copyright © 2009 by J. Stolfi, UNICAMP"
/* Created on 2009-01-17 by J.Stolfi, UNICAMP */
/* Last edited on 2024-12-05 10:40:24 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <dspmat.h>

void dspmat_linsys_SUBGE_solve
  ( double b[],
    dspmat_dim_t nb,
    dspmat_t *A, 
    double x[], 
    dspmat_dim_t nx,
    int32_t max_iter, 
    double omega,
    double abs_tol, 
    double rel_tol
  );
  /* Computes the solution {x} of the linear equation system {b = A x},
    where {A} is a given {nb × nx} sparse matrix, {b[0..nb-1]} is
    a given (column) vector, and {x[0..nx-1]} is an unknown (column)
    vector.
    
    Requires the matrix to be square, hence {nb == nx}. Said another
    way, stores into {x} the product {A^{-1}*b}.  The matrix must not
    contain two entries with the same indices.
    
    Uses a bootstrap method method with parameters
    {omega,abs_tol,rel_tol}, starting with {x = (0,0,.. 0)}. Stops
    after doing {max_iter} passes, or as soon as the residual is small
    according to {abs_tol,rel_tol}. */ 

void dspmat_inv_mul_SUBGE
  ( dspmat_t *A, 
    dspmat_t *B, 
    dspmat_t *X,
    int32_t max_iter,
    double abs_tol, 
    double rel_tol
  );
  /* Stores into {X} the product {A^{-1}*B}, namely, the solution of
    the linear equation system {A*X = B}.
    
    The matrices {A} and {B} must have dimensions {na × nx} and {na ×
    nb}, respectively, with {nx >= na}; the solution {X} will have
    dimensions {nx × nb}. Neither matrix may contain two entries with
    the same indices.
    
    Uses a partial elimination method to precondition the system,
    followed by Gauss-Seidel iteration starting with {x = (0,0,.. 0)}.
    Stops after doing {max_iter} passes, or as soon as the residual is
    small according to {abs_tol,rel_tol} */ 

void dspmat_linsys_SUBGE_precond(dspmat_t *A, dspmat_t *B);
  /* Modifies the matrices {A} and {B}, so as to turn the linear
    equation system {A*X = B} into an equivalent system for which the
    Gauss-Seidel iteration is guaranteed to converge.
    The requirements on the dimensions of {A} and {B}
    are those of {dspmat_linsys_SUBGE_solve}. */

#endif
