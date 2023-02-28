/* gauss_elim.h - Gaussian triangulation and elimination. */
/* Last edited on 2023-02-27 08:11:39 by stolfi */

#ifndef gauss_elim_H
#define gauss_elim_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

/* In all the procedures below, two-dimensional matrices are stored
  into one-dimensional vectors, in row-by-row order. That is, an {m×n}
  matrix {A[0..m-1,0..n-1]} is stored as a vector {A[0..m*n-1]}, 
  with entry {A[i,j]} of the matrix in element {A[n*i+j]} of the
  vector. */
    
int32_t gsel_solve(int32_t m, int32_t n, double A[], int32_t p, double B[], double X[], double tiny);
  /* Solves the linear system {A X = B}, where {A} is a known matrix
    of size {m × n}, {B} is a known matrix of size {m × p}, and {X} is
    an unknown matrix of size {n × p}. The arrays {A} and {B} are not
    modified.
    
    Returns the number {r} of equations used in the computation of
    {X}, in the range {0..min(m,n)}. If {r < m}, the matrix {A} has linearly
    dependent rows, and the result {X} is meaningless. If {r == m} but
    {m < n}, then the system has multiple solutions, and {n-m} rows of
    {X} have been set to zero.
    
    During triangularization, any entry whose absolute value gets
    reduced to {tiny} or less times its previous value is set to zero.
    If {tiny} is zero or negative, this cleanup is supressed. */

void gsel_residual(int32_t m, int32_t n, double A[], int32_t p, double B[], double X[], double R[]);
  /* Given the matrices {A} and {B} and a putative solution {X} of the system
    {A X = B}, computes the residual {R = A X - B}.  Assumes that {R}
    has size {m × p} (like {B}). */

double gsel_determinant(int32_t m, int32_t n, double A[], int32_t q);
  /* Returns the determinant of the first {q} rows and {q} columns of
    the {m × n} array {A}. If {q > m} or {q > n}, the result is zero.
    The array {A} is not changed. */

/* ARRAY TRIANGULARIZATION, DIAGONALIZATION, NORMALIZATION, DETERMINANT */

void gsel_triangularize(int32_t m, int32_t n, double M[], bool_t total, double tiny);
  /* Applies the Gaussian elimination method to the {m×n} matrix
    {M}, leaving it upper triangular. 
    
    Specifically, the matrix {M} is modified by row operations that do
    not change its determinant. Let {lead(M,i)} be the column index of
    the first non-zero element in row {i} of matrix {M}, or {+oo} if
    that row is all zeros. Upon exit, in any case, we will have
    {lead(M,i) >= i} for every {i}. Thus, in particular, if {m = n},
    the determinant of {M} will be the product of its diagonal
    elements.
    
    If {total} is true, the output matrix satisfies a stronger
    condition: for all {i < k}, {lead(M,i) < lead(M,k)}, or {lead(M,i)
    == lead(M,k) = +oo}. In this case, the non-zero rows of {M} are a
    basis for the space spanned by the original rows.
    
    During triangularization, any entry whose absolute value gets
    reduced to {tiny} or less times its previous value is set to zero.
    If {tiny} is zero or negative, this cleanup is supressed. */

void gsel_diagonalize(int32_t m, int32_t n, double M[]);
  /* Assumes that the {m×n} matrix {M} has been triangularized with
    {total = TRUE}, namely that {lead(M,i) = +oo} or {lead(M,i) >=
    lead(M,i-1)} for every {i}. Applies row operations to {M} so that it
    becomes diagonal-like, without change to its determinant. Namely,
    whenever {j := lead(M,i)} is finite, all elements {M[k,j]} in rows
    {k < i} are set to zero. */
    
void gsel_normalize(int32_t m, int32_t n, double M[]);
  /* Assumes that the {m×n} matrix {M} has been triangularized with
    {total=TRUE}, namely that {lead(M,i)} is {+oo} or {lead(M,i) >
    lead(M,i-1)} for every {i}; and possibly diagonalized).
    
    Then, whenever {j := lead(M,i)} is finite, scales the row
    so that {M[i,j] == 1.0}. */

/* SYSTEM SOLVING UTILITIES */

double gsel_triangular_det(int32_t m, int32_t n, double M[], int32_t q);
  /* Assumes that the {m×n} matrix {M} has been triangularized.
    Returns the product of the elements on the main diagonal of the
    first {q} rows and columns. Returns zero if {m < q} or {n < q}.
    Note that when {q == m} the result is the determinant of the first
    {m} columns of {M}. */

int32_t gsel_extract_solution(int32_t m, int32_t n, double M[], int32_t p, double X[]);
  /* Assumes that {M} is an {m×n} matrix that has been
    triangularized with {total = TRUE}, diagonalized, and normalized flag.
    
    The procedure interprets {M} as two blocks {M = (A B)} side by
    side, where {A} is {m×q} and {B} is {m×p}, with {q = n-p}. The
    procedure tries to compute a {q × p} matrix {X} such that {A X =
    B}. That is, it tries to solve the {p} linear systems {A x = u}
    where {x} is an unknown vector of size {q} and {u} is replaced in
    turn by each of the columns of {B}.
    
    Returns the number {r} of equations used in the computation of
    {X}, in the range {0..m}. If {r < m}, the matrix {Y = A X - B} may
    have {m-r} nonzero rows. If {r < q}, then the system has multiple
    solutions, and {q-m} rows of {X} have been set to zero.
    
    Specifically,
    
      * The returned result {r} is the number of rows {i < min(m,q)} such that
      {lead(M,i) < q}.
    
      * If {j} in {0..q-1} is not {lead(M,i)} for any
      {i}, then {X[j,0..p-1]} will be set to 0.
      
      * For every {i} in {0..m-1} such that {M[i,0..q-1] == 0},
      {Y[i,0..p-1]} will be equal to {-B[i,0..p-1]}.
    
    The matrix {X[i,k]} is stored in the given vector {X},
    linearized by rows. */

/* PRINTOUT */

void gsel_print_array(FILE *wr, char *fmt, char *head, int32_t m, int32_t n, double M[], char *foot);
  /* Writes to {wr} the array {M}, in a human-readable format. Each
    element is printed with format {fmt}. The strings {head} and
    {foot}, if not NULL, are printed on separate lines before and
    after the array. */

void gsel_print_system
  ( FILE *wr, 
    char *fmt, 
    char *head, 
    int32_t m, 
    int32_t n, 
    double A[], 
    int32_t p, 
    double B[], 
    char *foot
  );
  /* Writes to {wr} the coefficient matrix {A} and the right-hand-side
    matrix {B} of the linear equation system {A X = B}, in a
    human-readable format. Each element is printed with format {fmt}.
    The strings {head} and {foot}, if not NULL, are printed on
    separate lines before and after the system. */

#endif
