/* rmxn.h --- m by n matrices and operations on them */
/* Last edited on 2024-12-05 10:28:27 by stolfi */

#ifndef rmxn_H
#define rmxn_H

#include <stdio.h>
#include <stdint.h>
#include <rn.h>

/* Matrices are assumed to be stored linearized by rows, without gaps;
  so, in a matrix {A} with {m} rows and {n} columns, the logical 
  element {A[i,j]} is actually {A[n*i + j]}.
  
  Unless said otherwise, all output matrices and vectors must be 
  disjoint from each other from all the input ones. */

void rmxn_zero(uint32_t m, uint32_t n, double *M);
  /* Stores in {M} the null matrix. */

void rmxn_ident(uint32_t m, uint32_t n, double *M);
  /* Stores in {M} an identity matrix, truncated to {m} rows
    and {n}columns. */

void rmxn_copy(uint32_t m, uint32_t n, double *A, double *M);
  /* Copies the matrix {A} into {M}. Both matrices are assumed to have
    {m} rows and {n} columns. */

void rmxn_get_row(uint32_t m, uint32_t n, double *A, uint32_t i, double *r);
void rmxn_set_row(uint32_t m, uint32_t n, double *A, uint32_t i, double *r);
  /* These two procedures copy row {i} of matrix {A} to and from
    {r[0..n-1]}, respectively. The matrix {A} is assumed to have {m} rows
    and {n} columns. */

void rmxn_get_col(uint32_t m, uint32_t n, double *A, uint32_t j, double *r);
void rmxn_set_col(uint32_t m, uint32_t n, double *A, uint32_t j, double *r);
  /* These two procedures copy column {j} of matrix {A} to and from
    {r[0..m-1]}, respectively. The matrix {A} is assumed to have {m}
    rows and {n} columns. */

void rmxn_scale(uint32_t m, uint32_t n, double s, double *A, double *M);
  /* Stores in {M} the matrix {s*A}.  The matrix {M} may be the 
    same as {A}. */
    
void rmxn_add(uint32_t m, uint32_t n, double *A, double *B, double *M);
  /* Stores in {M} the matrix {A+B}.  The matrix {M} may be the 
    same as {A} and/or {B}. */
    
void rmxn_sub(uint32_t m, uint32_t n, double *A, double *B, double *M);
  /* Stores in {M} the matrix {A-B}.  The matrix {M} may be the 
    same as {A} and/or {B}. */

void rmxn_mix (uint32_t m, uint32_t n, double s, double *A, double t, double *B, double *M);
  /* Stores in {M} the matrix {s*A + t*B}. The matrix {M} may be the 
    same as {A}.*/

void rmxn_rel_diff(uint32_t m, uint32_t n, double *A, double *B, double *M);
  /* Stores in {z = M[i,j]} the relative difference between {x = A[i,j]}
    and {y = B[i,j]}, namely {(x-y)/sqrt(x^2 + y^2)/2}.  
    (See {rel_diff} in {jsmath.h}). The matrix {M} may be the 
    same as {A} or {B}. */

void rmxn_map_row (uint32_t m, uint32_t n, double *x, double *A, double *r);
  /* Computes {r = x * A}, where {x} is a (row) vector of size {m},
    and {A} is matrix of size {m x n}. The result is a (row) vector of
    size {n}.  The vector {r} must be disjoint from {A} and {x}. */

void rmxn_map_col (uint32_t m, uint32_t n, double *A, double *x, double *r);
  /* Computes {r = A * x}, where {A} is a matrix of size {m x n},
    and {x} is a (column) vector of size {n}. The result is a (column)
    vector of size {m}.  The vector {r} must be disjoint from {A} and {x}. */

void rmxn_mul (uint32_t m, uint32_t p, uint32_t n, double *A, double *B, double *M);
  /* Computes the matrix product {M = A * B}, where {A} has size
    {m x p} and {B} has size {p x n}. The matrix {M} must be disjoint
    from {A} and {B}. */

void rmxn_mul_tr (uint32_t m, uint32_t n, uint32_t p, double *A, double *B, double *M);
  /* Computes the matrix product {M = A * B^t}, where {A} has size
    {m x p} and {B} has size {n x p}.  (In other words, sets {M[i,j]}
    to the dot product of row {i} of {A} and row {j} of {B}.)
    The matrix {M} must have size {m x n} and must 
    be disjoint from {A} and {B}.  */

void rmxn_tr_mul (uint32_t p, uint32_t m, uint32_t n, double *A, double *B, double *M);
  /* Computes the matrix product {M = A^t * B}, where {A} has size
    {p x m} and {B} has size {p x n}.  (In other words, sets {M[i,j]}
    to the dot product of column {i} of {A} and column {j} of {B}.)
    The matrix {M} must have size {m x n} and must 
    be disjoint from {A} and {B}.  */

double rmxn_det (uint32_t n, double *A);
  /* Returns the determinant of the {n x n} matrix {A} */

#define rmxn_det_by_enum_SIZE_MAX 9
  /* Max value of {q} for {rmxn_det_by_enum}. Note that {9! = 362'880}. */

double rmxn_det_by_enum(uint32_t m, uint32_t n, double A[], uint32_t q);
  /* Determinant of the first {q} rows and columns of {A}, computed by
    the elementary definition (sum of {q!} products of elements of {A})
    Returns zero if {q > m} or {q > n}.  Otherwise {q}
    must not exceed {rmxn_det_by_enum_SIZE_MAX}. */

double rmxn_inv (uint32_t n, double *A, double *M);
  /* Sets {M} to the inverse of the {n x n} matrix {A}. The matrix {M}
    may be the same as {A}. Returns the determinant of {A}.
    
    If the returned determinant is zero, the contents of {M} should be
    considered garbage. Note that if {A} is singular or nearly so the
    procedure may still return a nonzero value, but the contents of
    {M} will probably be garbage all the same. 
    
    The current implementation uses Gaussian elimination with full
    pivoting. */

double rmxn_norm_sqr(uint32_t m, uint32_t n, double *A);
  /* Squared Frobenius norm of the {m � n} matrix {A}, 
    i.e. sum of squares of elements */

double rmxn_norm(uint32_t m, uint32_t n, double *A);
  /* Frobenius norm of the {m � n} matrix {A},
    i.e. square root of sum of squares of elements */

double rmxn_normalize(uint32_t m, uint32_t n, double *A);
  /* Divides the matrix {A} by its {rmxn_norm}. 
    If that norm is zero, the matrix is filled with {NAN}.
    Returns the norm. */

double rmxn_mod_norm_sqr (uint32_t n, double *A);
  /* Returns the square of the Frobenius norm of {A-I}.
    The matrix {A} is assumed to be {n � n}. */

double rmxn_max_abs_elem(uint32_t m, uint32_t n, double *A);
  /* Returns the maximum absolute value of any element in the 
    {m � n} matrix {A}. */

double rmxn_max_abs_elem_in_row(uint32_t m, uint32_t n, double M[], uint32_t i);
  /* Returns the maximum absolute value of the elements in row {i}
    of the {m � n} matrix {M}. */

double rmxn_max_abs_elem_in_col(uint32_t m, uint32_t n, double M[], uint32_t j);
  /* Returns the maximum absolute value of the elements in column {j} of
     the {m � n} matrix {M}. */
    
/* MATRIX FACTORIZATION */

void rmxn_cholesky(uint32_t n, double *A, double *L);
  /* Factors the positive definite matrix {A}, with {n} rows and
   columns, into {L*L^t} where {L} is lower triangular and {L^t} is
   the transpose of {L}. The matrix {L} must be distinct from {A}.
   
   (The name is French, thus it should be pronounced "sholesKEE"). */

/* OPERATIONS ON LOWER TRIANGULAR MATRICES */

/* The operations in this section require the matrix {L} to 
   be lower triangular, with nonzero entries in the diagonal. */

void rmxn_LT_inv_map_row(uint32_t n, double *y, double *L, double *r);
  /* Computes {r = y * (L^-1)}, where {y} is a (row) vector of size
    {n}, and {L} is a lower triangular matrix of size {n} by {n}. The
    result is a (row) vector of size {n}. In other words, finds the
    solution {r} to the linear system {r*L = y}. The vector {r} may be
    the same as {y}. */
    
void rmxn_LT_inv_map_col(uint32_t m, double *L, double *y, double *r);
  /* Computes {r = (L^-1) * y}, where {L} is a lower triangular matrix
    of size {m} by {m}, and {y} is a (column) vector of size {m}. The
    result is a (column) vector of size {m}. In other words, finds the
    solution {r} to the linear system {L*r = y}. The vector {r} may be
    the same as {y}. */

void rmxn_LT_pos_div(uint32_t m, uint32_t n, double *A, double *L, double *M);
  /* Computes {M = A * (L^-1)}, where {A} is a rectangular matrix of
    size {m} by {n}, and {L} is a lower triangular matrix of size {n}
    by {n}. In other words, finds the solution {M} to the matrix
    equation {M*L = A}. The result {M} has size {m} by {n}, and may be
    the same as {A}. */

void rmxn_LT_pre_div(uint32_t m, uint32_t n, double *L, double *A, double *M);
  /* Computes {M = (L^-1) * A}, where {L} is a lower triangular matrix
    of size {m} by {m}, and {A} is a rectangular matrix of size
    {m} by {n}. In other words, finds the solution {M} to the matrix
    equation {L*M = A}. The result has size {m} by {n}, and may be
    the same as {A}. */

void rmxn_perturb_unif(uint32_t m, uint32_t n, double pabs, double prel, double *A);
  /* Adds an independent random perturbation with uniform
    distribution in {[-mag _ +mag]} to each element {A[i][j]} of {M},
    where {mag} is {pabs + prel*fabs(M[i][j])}. The matrix {A} is
    assumed to have {m*n} elements. */ 
    
void rmxn_cleanup(uint32_t m, uint32_t n, double *A, double tiny);
  /* Sets to zero any elements of {A} that is less than {tiny} in 
    absolute value.  A no-op if {tiny} is zero or negative. */

/* SIMPLE MATRIX PRINTOUT */

void rmxn_print (FILE *f, uint32_t m, uint32_t n, double *A);
  /* Prints the {m x n} matrix {A} to file {f}, with default format. */
    
/* FLEXIBLE MATRIX PRINTOUT 

  The following procedures use {fmt} to print each element. Each matrix
  is bounded by {olp} and {orp}, and rows are separated by {osep}. Each
  row is bounded by {ilp} and {irp}, and elements are separated by
  {isep}. If there are two or more matrices side by side, the string
  {msep} is printed between them, in each row. Defaults are provided for
  any of {olp,osep,orp,ilp,isep,irp,msep} which are NULL.
  
  If {m} is zero, only {olp,osep,orp} are printed, and no matrix rows
  are printed.
  
  If {m} is positive but the number of columns of any matrix
  ({n},{n1},{n2}, or {n3}) is zero, an empty matrix of {m} rows (with
  the {ilp} and {irp} delimiters only) is printed. If any of the
  matrices ({A},{A1},{A2},{A3}) is NULL, the matrix and its preceding
  {msep}, if any, are omitted. */

void rmxn_gen_print 
  ( FILE *f, 
    uint32_t m, 
    uint32_t n, double *A,
    char *fmt, 
    char *olp, char *osep, char *orp, /* Outer delimiters. */
    char *ilp, char *isep, char *irp  /* Inner delimiters. */
  );
  /* Prints the {m x n} matrix {A} to file {f}. 
    The */

void rmxn_gen_print2 
  ( FILE *f, 
    uint32_t m,
    uint32_t n1, double *A1,
    uint32_t n2, double *A2,
    char *fmt, 
    char *olp, char *osep, char *orp,     /* Outer delimiters. */
    char *ilp, char *isep, char *irp,     /* Inner delimiters. */
    char *msep                            /* Matrix separator. */
  );
  /* Prints the {m x n1} matrix {A1} and the {m x n2} matrix {A2}, side by side, to file {f}. */
 
void rmxn_gen_print3 
  ( FILE *f, 
    uint32_t m,
    uint32_t n1, double *A1,
    uint32_t n2, double *A2,
    uint32_t n3, double *A3,
    char *fmt, 
    char *olp, char *osep, char *orp,     /* Outer delimiters. */
    char *ilp, char *isep, char *irp,     /* Inner delimiters. */
    char *msep                            /* Matrix separator. */
  );
  /* Prints the {m x n1} matrix {A1}, the {m x n2} matrix {A2},
    and the {m x n3} matrix {A3}, side by side, to file {f}. */

/* HEAP ALLOCATION */

double *rmxn_alloc(uint32_t m, uint32_t n);
  /* Allocates {m*n} {double}s on the heap; bombs out if no mem. */

#endif

