/* rmxn_throw.h --- random MxN matrices */
/* Last edited on 2024-12-05 10:28:44 by stolfi */

#ifndef rmxn_throw_H
#define rmxn_throw_H

#include <stdint.h>

void rmxn_throw_matrix(uint32_t m, uint32_t n, double *A);
  /* Fills {A} with an {m} by {n} matrix whose elements
    are all independent random numbers in the range {[-1 _ +1]}. */

void rmxn_throw_singular(uint32_t n, double *A);
  /* Fills {A} with a random {n} matrix with unit RMS element value
    whose determinant is zero apart form roundoff errors. 
    The size {n} must be at least 2. */
    
void rmxn_throw_almost_singular_pair(uint32_t n, double *A, double *B, double detMax);
  /* Fills the {n} by {n} matrices {A} and {B} with random numbers such
    that both {A} and {B} have unit RMS element size and nonzero
    determinants, and at least one of the determinants has absolute
    value at most {detMax}, and the product {A*B} is the identity apart
    from a positive homogeneous scaling factor and roundoff errors.
    
    In particular, if {n==1} then {detMax} must be at least 1.0, and 
    both {A and {B} are set to the {1x1} ideitity matrix.
    
    The procedure may take a long time, or fail altogether, if the
    parameter {detMax} is too small. */
    
void rmxn_throw_non_singular_pair(uint32_t n, double *A, double *B, double detMin);
  /* Fills the {n} by {n} matrices {A} and {B} with random numbers such
    that both {A} and {B} have unit RMS element size, both have
    determinants with absolute value at least {detMin}, and the product
    {A*B} is the identity apart from a positive homogeneous scaling
    factor and roundoff errors.
    
    In particular, if {n==1} then {detMin} must be at most 1.0, and 
    both {A and {B} are set to the {1x1} ideitity matrix.
    
    The procedure may take a long time, or fail altogether, if the
    parameter {detMin} is too large. */

void rmxn_throw_ortho(uint32_t n, double *A);
  /* Stores in {A} a random positive orthonormal {n×n} matrix. 
    Its rows and columns are pairwise orthogonal and with length 1, and ditto
    for its columns; and its determinant is {+1}.  Assumes {A} has {n^2} elements. */

void rmxn_throw_ortho_complement(uint32_t n, uint32_t ma, double *A, uint32_t mb, double *B);
  /* Assumes that {A} is a {ma} by {n} matrix and {B} is a {mb} by {n}
    matrix, both stored by rows, with {ma + d <= n}. Assumes that the
    rows of {A} are orthonormal. Stores into the rows of {B} a set of
    {mb} {n}-vectors that are orthonormal and orthogonal to all rows of
    {A}. If {ma} is zero then {A} may be {NULL}. */

void rmxn_throw_directions(uint32_t m, uint32_t n, double U[]);
  /* Assumes that {U} is an {m} by {n} matrix, stored by rows.
    Stores into the rows of {U} a set of {m} well-spaced directions
    in {\RR^n}.  
    
    For {i} in {0..min(m,n)}, row {i} is the direction
    of coordinate axis {i}.  If {m > n}, the following rows are
    random directions that make a sufficiently large angle 
    with the previous rows.
    
    If {n} is 0 or 1, {m} must be equal to {n}. */

void rmxn_throw_LT_matrix(uint32_t m, double *Lmm);
  /* Fills {Lmm} with a lower triangular matrix where element {Lmm[i,j]}
    is a random value in {[-1 _ +1]} if {j <= i}, and zero otherwise. */

#endif

