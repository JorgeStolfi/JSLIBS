#ifndef dspmat_linsys_GS_H
#define dspmat_linsys_GS_H
/* Gauss-Seidel linear system solving for sparse matrices with {double} entries */

#define dspmat_linsys_GS_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-08-17 by J.Stolfi, UNICAMP */
/* Last edited on 2023-03-18 10:49:07 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <dspmat.h>

void dspmat_linsys_GS_solve
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
    
    Uses the iterative Gauss-Seidel method (see
    {dspmat_linsys_GS_solve_iteration}) with parameters
    {omega,abs_tol,rel_tol}, starting with {x = (0,0,.. 0)}. Stops
    after doing {max_iter} passes, or as soon as
    {dspmat_linsys_GS_solve_iteration} returns a result less than 1.0. */ 

void dspmat_inv_mul_GS
  ( dspmat_t *A, 
    dspmat_t *B, 
    dspmat_t *X,
    int32_t max_iter, 
    double omega,
    double abs_tol, 
    double rel_tol
  );
  /* Sets {*X} to {A^{-1}*B}, the inverse of {A} times the matrix {B}; namely,
    the solution of the matrix equation {A*X == B}. The matrix {A} must be
    square.
    
    Uses {dspmat_linsys_GS_solve}) with parameters
    {max_iter,omega,abs_tol,rel_tol}. Basically equivalent to solving
    {b = A*x} where {b} is each column of {B} and {x} is the
    corresponding column of {X}. */ 

double dspmat_linsys_GS_solve_iteration
  ( double b[],
    dspmat_t *A, 
    double x[],
    double omega,
    double abs_tol, 
    double rel_tol
  );
  /* 
    Given an {n × n} matrix {A} and two {n}-vectors {b,x}, performs
    one iteration of the Gauss-Seidel algorithm for the linear
    equation system {b=A*x}, with relaxation factor {omega}. Namely,
    for each {i} in {0..n-1}, computes the quantity
    
      { e[i] = b[i] - SUM{A[i,j]*x[j] : j\neq i}}
      
    then computes the new value { xNew[i] = e[i]/A[i,i] },
    then sets {x[i] = omega*xNew[i] + (1-omega)*x[i]}.
    
    Returns the maximum residual {b[i] = e[i] - A[i,i]*x[i]} observed
    during this computation, divided by a `tolerance criterion' which
    is a combination of {abs_tol}, AND {rel_tol} times the the
    magnitude of {b[i]}. Thus, if the result is less than 1, one may
    assume that the residual was `small enough'. */

#endif
