/* rmxn_extra.h --- additional operation on MxN matrices */
/* Last edited on 2024-11-07 22:36:23 by stolfi */

#ifndef rmxn_extra_H
#define rmxn_extra_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

void rmxn_throw(int32_t m, int32_t n, double M[]);
  /* Fills {M} with random numbers in the interval {[-1 _ +1]}. */

void rmxn_perturb_unif(int32_t m, int32_t n, double pabs, double prel, double M[]);
  /* Adds an independent random perturbation of with uniform distribution 
    in {[-mag _ +mag]} to each element {M[i][j]} of {M}, where {mag} is {pabs + prel*fabs(M[i][j])}.
    The matrix {M} is assumed to have {m*n} elements. */ 

void rmxn_throw_ortho(int32_t n, double M[]);
  /* Stores in {M} a random positive orthonormal {n×n} matrix. 
    Its rows and columns are pairwise orthogonal and with length 1, and ditto
    for its columns; and its determinant is {+1}.  Assumes {M} has {n^2} elements. */

void rmxn_throw_ortho_complement(int32_t n, int32_t p, double A[], int32_t q, double M[]);
  /* Assumes that {A} is a {p} by {n} matrix and {M} is a {q} by {n} matrix, 
    both stored by rows, with {p + d <= n}. Assumes that the rows of {A} are orthonormal.
    Stores into the rows of {M} a set of {q} {n}-vectors that are orthonormal and orthogonal
    to all rows of {A}. If {p} is zero then {A} may be {NULL}. */

void rmxn_throw_directions(int32_t m, int32_t n, double U[]);
  /* Assumes that {U} is an {m} by {n} matrix, stored by rows.
    Stores into the rows of {U} a set of {m} well-spaced directions
    in {\RR^n}.  
    
    For {i} in {0..min(m,n)}, row {i} is the direction
    of coordinate axis {i}.  If {m > n}, the following rows are
    random directions that make a sufficiently large angle 
    with the previous rows. */

void rmxn_spin_rows(int32_t m, int32_t n, double A[], double M[]);
  /* Applies a random rotation to each row of {A}, which is assumed to
    have {m} rows and {n} columns. Equivalent to computing {M=A*N}
    where {N} is a random orthonormal {n×n} matrix. */

void rmxn_spin_cols(int32_t m, int32_t n, double A[], double M[]);
  /* Applies a random rotation to each column {A}, which is assumed to
    have {m} rows and {n} columns. Equivalent to computing {M=N*A}
    where {N} is a random orthonormal {m×m} matrix. */

void rmxn_shift_rows(int32_t m, int32_t n, double A[], double v[], double M[]);
  /* Adds the vector {v[0..n-1]} to each row of the matrix {A},
    yielding matrix {M}. Both matrices are assumed to have {m} rows
    and {n} columns. */

void rmxn_shift_cols(int32_t m, int32_t n, double A[], double v[], double M[]);
  /* Adds the vector {v[0..m-1]} to each column of the matrix {A},
    yielding matrix {M}. Both matrices are assumed to have {m} rows
    and {n} columns. */

#endif

