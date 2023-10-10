/* rmxn.h --- m by n matrices and operations on them */
/* Last edited on 2023-10-09 09:03:31 by stolfi */

#ifndef rmxn_H
#define rmxn_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <rn.h>

/* Matrices are assumed to be stored linearized by rows, without gaps;
  so, in a matrix {A} with {m} rows and {n} columns, the logical 
  element {A[i,j]} is actually {A[n*i + j]}.
  
  Unless said otherwise, all output matrices and vectors must be 
  disjoint from each other from all the input ones. */

void rmxn_zero(int32_t m, int32_t n, double *M);
  /* Stores in {M} the null matrix. */

void rmxn_ident(int32_t m, int32_t n, double *M);
  /* Stores in {M} an identity matrix, truncated to {m} rows
    and {n}columns. */

void rmxn_copy(int32_t m, int32_t n, double *A, double *M);
  /* Copies the matrix {A} into {M}. Both matrices are assumed to have
    {m} rows and {n} columns. */

void rmxn_get_row(int32_t m, int32_t n, double *A, int32_t i, double *r);
void rmxn_set_row(int32_t m, int32_t n, double *A, int32_t i, double *r);
  /* These two procedures copy row {i} of matrix {A} to and from
    {r[0..n-1]}, respectively. The matrix {A} is assumed to have {m} rows
    and {n} columns. */

void rmxn_get_col(int32_t m, int32_t n, double *A, int32_t j, double *r);
void rmxn_set_col(int32_t m, int32_t n, double *A, int32_t j, double *r);
  /* These two procedures copy column {j} of matrix {A} to and from
    {r[0..m-1]}, respectively. The matrix {A} is assumed to have {m}
    rows and {n} columns. */

void rmxn_scale(int32_t m, int32_t n, double s, double *A, double *M);
  /* Stores in {M} the matrix {s*A}.  The matrix {M} may be the 
    same as {A}. */

void rmxn_mix (int32_t m, int32_t n, double s, double *A, double t, double *B, double *M);
  /* Stores in {M} the matrix {s*A + t*B}. The matrix {M} may be the 
    same as {A}.*/

void rmxn_rel_diff(int32_t m, int32_t n, double *A, double *B, double *M);
  /* Stores in {z = M[i,j]} the relative difference between {x = A[i,j]}
    and {y = B[i,j]}, namely {(x-y)/sqrt(x^2 + y^2)/2}.  
    (See {rel_diff} in {jsmath.h}). The matrix {M} may be the 
    same as {A} or {B}. */

void rmxn_map_row (int32_t m, int32_t n, double *x, double *A, double *r);
  /* Computes {r = x * A}, where {x} is a (row) vector of size {m},
    and {A} is matrix of size {m x n}. The result is a (row) vector of
    size {n}.  The vector {r} must be disjoint from {A} and {x}. */

void rmxn_map_col (int32_t m, int32_t n, double *A, double *x, double *r);
  /* Computes {r = A * x}, where {A} is a matrix of size {m x n},
    and {x} is a (column) vector of size {n}. The result is a (column)
    vector of size {m}.  The vector {r} must be disjoint from {A} and {x}. */

void rmxn_mul (int32_t m, int32_t p, int32_t n, double *A, double *B, double *M);
  /* Computes the matrix product {M = A * B}, where {A} has size
    {m x p} and {B} has size {p x n}. The matrix {M} must be disjoint
    from {A} and {B}. */

void rmxn_mul_tr (int32_t m, int32_t n, int32_t p, double *A, double *B, double *M);
  /* Computes the matrix product {M = A * B^t}, where {A} has size
    {m x p} and {B} has size {n x p}.  (In other words, sets {M[i,j]}
    to the dot product of row {i} of {A} and row {j} of {B}.)
    The matrix {M} must have size {m x n} and must 
    be disjoint from {A} and {B}.  */

void rmxn_tr_mul (int32_t p, int32_t m, int32_t n, double *A, double *B, double *M);
  /* Computes the matrix product {M = A^t * B}, where {A} has size
    {p x m} and {B} has size {p x n}.  (In other words, sets {M[i,j]}
    to the dot product of column {i} of {A} and column {j} of {B}.)
    The matrix {M} must have size {m x n} and must 
    be disjoint from {A} and {B}.  */

double rmxn_det (int32_t n, double *A);
  /* Returns the determinant of the {n x n} matrix {A} */

double rmxn_cof (int32_t n, double *A, int32_t ix, int32_t jx);
  /* Returns the cofactor of element {[ix,jx]} in the {n x n} matrix {A}. */

void rmxn_adj (int32_t n, double *A, double *M);
  /* Sets {M} to the adjoint of the {n x n} matrix {A}.
    The matrix {M} may be the same as {A}. */

double rmxn_inv (int32_t n, double *A, double *M);
double rmxn_inv_full (int32_t n, double *A, double *M);
  /* Sets {M} to the inverse of the {n x n} matrix {A}. The matrix {M}
    may be the same as {A}. Returns the determinant of {A}.
    
    If the returned determinant is zero, the contents of {M} should be
    considered garbage. Note that if {A} is singular or nearly so the
    procedure may still return a nonzero value, but the contents of
    {M} will probably be garbage all the same. The version
    {rmxn_inf_full} uses full pivoting and therefore may be more
    robust and/or accurate in these cases. */

double rmxn_norm_sqr(int32_t m, int32_t n, double *A);
  /* Squared Frobenius norm of the {m × n} matrix {A}, 
    i.e. sum of squares of elements */

double rmxn_norm(int32_t m, int32_t n, double *A);
  /* Frobenius norm of the {m × n} matrix {A},
    i.e. square root of sum of squares of elements */

double rmxn_normalize(int32_t m, int32_t n, double *A);
  /* Divides the matrix {A} by its {rmxn_norm}. 
    If that norm is zero, the matrix is filled with {NAN}.
    Returns the norm. */

double rmxn_mod_norm_sqr (int32_t n, double *A);
  /* Returns the square of the Frobenius norm of {A-I}.
    The matrix {A} is assumed to be {n × n}. */

double rmxn_max_abs_elem(int32_t m, int32_t n, double *A);
  /* Returns the maximum absolute value of any element in the 
    {m × n} matrix {A}. */
    
/* MATRIX FACTORIZATION */

void rmxn_cholesky(int32_t n, double *A, double *L);
  /* Factors the positive definite matrix {A}, with {n} rows and
   columns, into {L*L^t} where {L} is lower triangular and {L^t} is
   the transpose of {L}. The matrix {L} must be distinct from {A}.
   
   (The name is French, thus it should be pronounced "sholesKEE"). */

/* OPERATIONS ON LOWER TRIANGULAR MATRICES */

/* The operations in this section require the matrix {L} to 
   be lower triangular, with nonzero entries in the diagonal. */

void rmxn_LT_inv_map_row(int32_t n, double *y, double *L, double *r);
  /* Computes {r = y * (L^-1)}, where {y} is a (row) vector of size
    {n}, and {L} is a lower triangular matrix of size {n} by {n}. The
    result is a (row) vector of size {n}. In other words, finds the
    solution {r} to the linear system {r*L = y}. The vector {r} may be
    the same as {y}. */
    
void rmxn_LT_inv_map_col(int32_t m, double *L, double *y, double *r);
  /* Computes {r = (L^-1) * y}, where {L} is a lower triangular matrix
    of size {m} by {m}, and {y} is a (column) vector of size {m}. The
    result is a (column) vector of size {m}. In other words, finds the
    solution {r} to the linear system {L*r = y}. The vector {r} may be
    the same as {y}. */

void rmxn_LT_pos_div(int32_t m, int32_t n, double *A, double *L, double *M);
  /* Computes {M = A * (L^-1)}, where {A} is a rectangular matrix of
    size {m} by {n}, and {L} is a lower triangular matrix of size {n}
    by {n}. In other words, finds the solution {M} to the matrix
    equation {M*L = A}. The result {M} has size {m} by {n}, and may be
    the same as {A}. */

void rmxn_LT_pre_div(int32_t m, int32_t n, double *L, double *A, double *M);
  /* Computes {M = (L^-1) * A}, where {L} is a lower triangular matrix
    of size {m} by {m}, and {A} is a rectangular matrix of size
    {m} by {n}. In other words, finds the solution {M} to the matrix
    equation {L*M = A}. The result has size {m} by {n}, and may be
    the same as {A}. */
    
/* SIMPLE MATRIX PRINTOUT */

void rmxn_print (FILE *f, int32_t m, int32_t n, double *A);
  /* Prints the {m x n} matrix {A} to file {f}, with default format. */
    
/* FLEXIBLE MATRIX PRINTOUT 

  The following procedures use {fmt} to print each element. Each matrix
  is bounded by {olp} and {orp}, and rows are separated by {osep}. Each
  row is bounded by {ilp} and {irp}, and elements are separated by
  {isep}. If there are two or more matrices side by side, the string
  {msep} is printed between them, in each row. Defaults are provided for
  any of {olp,osep,orp,ilp,isep,irp,msep} which are NULL. 
  
  If the number of columns ({n},{n1},{n2}, or {n3}) is given as 0, the
  corresponding matrix may be {NULL}, and an empty matrix (with the
  {ilp} and {irp} delimiters only) is printed. If the number of columns
  is negative, the matrix and its preceding {msep}, if any, are
  omitted. */

void rmxn_gen_print 
  ( FILE *f, int32_t m, int32_t n, double *A,
    char *fmt, 
    char *olp, char *osep, char *orp, /* Outer delimiters. */
    char *ilp, char *isep, char *irp  /* Inner delimiters. */
  );
  /* Prints the {m x n} matrix {A} to file {f}. */

void rmxn_gen_print2 
  ( FILE *f, int32_t m,
    int32_t n1, double *A1,
    int32_t n2, double *A2,
    char *fmt, 
    char *olp, char *osep, char *orp,     /* Outer delimiters. */
    char *ilp, char *isep, char *irp,     /* Inner delimiters. */
    char *msep                            /* Matrix separator. */
  );
  /* Prints the {m x n1} matrix {A1} and the {m x n2} matrix {A2}, side by side, to file {f}. */
 
void rmxn_gen_print3 
  ( FILE *f, int32_t m,
    int32_t n1, double *A1,
    int32_t n2, double *A2,
    int32_t n3, double *A3,
    char *fmt, 
    char *olp, char *osep, char *orp,     /* Outer delimiters. */
    char *ilp, char *isep, char *irp,     /* Inner delimiters. */
    char *msep                            /* Matrix separator. */
  );
  /* Prints the {m x n1} matrix {A1}, the {m x n2} matrix {A2}, and the {m x n3} matrix {A3}, side by side, to file {f}. */

/* HEAP ALLOCATION */

double *rmxn_alloc(int32_t m, int32_t n);
  /* Allocates {m*n} {double}s on the heap; bombs out if no mem. */

#endif

