/* rfmxn.h --- m by n matrices and operations on them */
/* Last edited on 2021-08-17 05:55:27 by stolfi */

#ifndef rfmxn_H
#define rfmxn_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <rn.h>

/* Matrices are assumed to be stored linearized by rows, without gaps;
  so, in a matrix {A} with {m} rows and {n} columns, the logical 
  element {A[i,j]} is actually {A[n*i + j]}.
  
  Unless said otherwise, all output matrices and vectors must be 
  disjoint from each other from all the input ones. */

void rfmxn_zero(int32_t m, int32_t n, float *M);
  /* Stores in {M} the null matrix. */

void rfmxn_ident(int32_t m, int32_t n, float *M);
  /* Stores in {M} an identity matrix, truncated to {m} rows
    and {n}columns. */

void rfmxn_copy(int32_t m, int32_t n, float *A, float *M);
  /* Copies the matrix {A} into {M}. Both matrices are assumed to have
    {m} rows and {n} columns. */

void rfmxn_get_row(int32_t m, int32_t n, float *A, int32_t i, float *r);
void rfmxn_set_row(int32_t m, int32_t n, float *A, int32_t i, float *r);
  /* These two procedures copy row {i} of matrix {A} to and from
    {r[0..n-1]}, respectively. The matrix {A} is assumed to have {m} rows
    and {n} columns. */

void rfmxn_get_col(int32_t m, int32_t n, float *A, int32_t j, float *r);
void rfmxn_set_col(int32_t m, int32_t n, float *A, int32_t j, float *r);
  /* These two procedures copy column {j} of matrix {A} to and from
    {r[0..m-1]}, respectively. The matrix {A} is assumed to have {m}
    rows and {n} columns. */

void rfmxn_scale(int32_t m, int32_t n, float s, float *A, float *M);
  /* Stores in {M} the matrix {s*A}.  The matrix {M} may be the 
    same as {A}. */

void rfmxn_mix (int32_t m, int32_t n, double s, float *A, double t, float *B, float *M);
  /* Stores in {M} the matrix {s*A + t*B}. The matrix {M} may be the 
    same as {A}.*/

void rfmxn_rel_diff(int32_t m, int32_t n, float *A, float *B, float *M);
  /* Stores in {z = M[i,j]} the relative difference between {x = A[i,j]}
    and {y = B[i,j]}, namely {(x-y)/sqrt(x^2 + y^2)/2}.  
    (See {rel_diff} in {jsmath.h}). The matrix {M} may be the 
    same as {A} or {B}. */

void rfmxn_map_row (int32_t m, int32_t n, float *x, float *A, float *r);
  /* Computes {r = x * A}, where {x} is a (row) vector of size {m},
    and {A} is matrix of size {m x n}. The result is a (row) vector of
    size {n}.  The vector {r} must be disjoint from {A} and {x}. */

void rfmxn_map_col (int32_t m, int32_t n, float *A, float *x, float *r);
  /* Computes {r = A * x}, where {A} is a matrix of size {m x n},
    and {x} is a (column) vector of size {n}. The result is a (column)
    vector of size {m}.  The vector {r} must be disjoint from {A} and {x}. */

void rfmxn_mul (int32_t m, int32_t p, int32_t n, float *A, float *B, float *M);
  /* Computes the matrix product {M = A * B}, where {A} has size
    {m x p} and {B} has size {p x n}. The matrix {M} must be disjoint
    from {A} and {B}. */

void rfmxn_mul_tr (int32_t m, int32_t n, int32_t p, float *A, float *B, float *M);
  /* Computes the matrix product {M = A * B^t}, where {A} has size
    {m x p} and {B} has size {n x p}.  (In other words, sets {M[i,j]}
    to the dot product of row {i} of {A} and row {j} of {B}.)
    The matrix {M} must have size {m x n} and must 
    be disjoint from {A} and {B}.  */

void rfmxn_tr_mul (int32_t p, int32_t m, int32_t n, float *A, float *B, float *M);
  /* Computes the matrix product {M = A^t * B}, where {A} has size
    {p x m} and {B} has size {p x n}.  (In other words, sets {M[i,j]}
    to the dot product of column {i} of {A} and column {j} of {B}.)
    The matrix {M} must have size {m x n} and must 
    be disjoint from {A} and {B}.  */

double rfmxn_det (int32_t n, float *A);
  /* Returns the determinant of the {n x n} matrix {A} */

double rfmxn_cof (int32_t n, float *A, int32_t ix, int32_t jx);
  /* Returns the cofactor of element {[ix,jx]} in the {n x n} matrix {A}. */

void rfmxn_adj (int32_t n, float *A, float *M);
  /* Sets {M} to the adjoint of the {n x n} matrix {A}.
    The matrix {M} may be the same as {A}. */

double rfmxn_inv (int32_t n, float *A, float *M);
double rfmxn_inv_full (int32_t n, float *A, float *M);
  /* Sets {M} to the inverse of the {n x n} matrix {A}. The matrix {M}
    may be the same as {A}. Returns the determinant of {A}.
    
    If the returned determinant is zero, the contents of {M} should be
    considered garbage. Note that if {A} is singular or nearly so the
    procedure may still return a nonzero value, but the contents of
    {M} will probably be garbage all the same. The version
    {rfmxn_inf_full} uses full pivoting and therefore may be more
    robust and/or accurate in these cases. */

double rfmxn_norm_sqr(int32_t m, int32_t n, float *A);
  /* Squared Frobenius norm of the {m × n} matrix {A}, 
    i.e. sum of squares of elements */

double rfmxn_norm(int32_t m, int32_t n, float *A);
  /* Frobenius norm of the {m × n} matrix {A},
    i.e. square root of sum of squares of elements */

double rfmxn_mod_norm_sqr (int32_t n, float *A);
  /* Returns the square of the Frobenius norm of {A-I}.
    The matrix {A} is assumed to be {n × n}. */

double rfmxn_max_abs_elem(int32_t m, int32_t n, float *A);
  /* Returns the maximum absolute value of any element in the 
    {m × n} matrix {A}. */
    
/* MATRIX FACTORIZATION */

void rfmxn_cholesky(int32_t n, float *A, float *L);
  /* Factors the positive definite matrix {A}, with {n} rows and
   columns, into {L*L^t} where {L} is lower triangular and {L^t} is
   the transpose of {L}. The matrix {L} must be distinct from {A}.
   
   (The name is French, thus it should be pronounced "sholesKEE"). */

/* OPERATIONS ON LOWER TRIANGULAR MATRICES */

/* The operations in this section require the matrix {L} to 
   be lower triangular, with nonzero entries in the diagonal. */

void rfmxn_LT_inv_map_row(int32_t n, float *y, float *L, float *r);
  /* Computes {r = y * (L^-1)}, where {y} is a (row) vector of size
    {n}, and {L} is a lower triangular matrix of size {n} by {n}. The
    result is a (row) vector of size {n}. In other words, finds the
    solution {r} to the linear system {r*L = y}. The vector {r} may be
    the same as {y}. */
    
void rfmxn_LT_inv_map_col(int32_t m, float *L, float *y, float *r);
  /* Computes {r = (L^-1) * y}, where {L} is a lower triangular matrix
    of size {m} by {m}, and {y} is a (column) vector of size {m}. The
    result is a (column) vector of size {m}. In other words, finds the
    solution {r} to the linear system {L*r = y}. The vector {r} may be
    the same as {y}. */

void rfmxn_LT_pos_div(int32_t m, int32_t n, float *A, float *L, float *M);
  /* Computes {M = A * (L^-1)}, where {A} is a rectangular matrix of
    size {m} by {n}, and {L} is a lower triangular matrix of size {n}
    by {n}. In other words, finds the solution {M} to the matrix
    equation {M*L = A}. The result {M} has size {m} by {n}, and may be
    the same as {A}. */

void rfmxn_LT_pre_div(int32_t m, int32_t n, float *L, float *A, float *M);
  /* Computes {M = (L^-1) * A}, where {L} is a lower triangular matrix
    of size {m} by {m}, and {A} is a rectangular matrix of size
    {m} by {n}. In other words, finds the solution {M} to the matrix
    equation {L*M = A}. The result has size {m} by {n}, and may be
    the same as {A}. */
    
/* MATRIX FORMATTING */

void rfmxn_print (FILE *f, int32_t m, int32_t n, float *A);
  /* Prints the {m x n} matrix {A} to file {f}, with default format. */

void rfmxn_gen_print 
  ( FILE *f, int32_t m, int32_t n, float *A,
    char *fmt, 
    char *olp, char *osep, char *orp, /* Outer delimiters. */
    char *ilp, char *isep, char *irp  /* Inner delimiters. */
  );
  /* Prints the {m x n} matrix {A} to file {f}, using {fmt} for each
    element. The matrix is bounded by {olp} and {orp}, and rows are
    separated by {osep}. Each row is bounded by {ilp} and {irp}, and
    elements are separated by {isep}. Defaults are provided for any of
    these strings which are NULL. */

/* HEAP ALLOCATION */

float *rfmxn_alloc(int32_t m, int32_t n);
  /* Allocates {m*n} {float}s on the heap; bombs out if no mem. */

#endif

