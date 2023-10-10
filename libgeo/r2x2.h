/* r2x2.h --- 2x2 matrices and operations on them */
/* Last edited on 2023-10-09 09:00:03 by stolfi */

#ifndef r2x2_H
#define r2x2_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <r2.h>

typedef struct r2x2_t { double c[2][2]; } r2x2_t;

/* STANDARD OPERATIONS */

void r2x2_zero(r2x2_t *M);
  /* Stores in {M} the null matrix. */

void r2x2_ident(r2x2_t *M);
  /* Stores in {M} the identity matrix. */

void r2x2_transp (r2x2_t *A, r2x2_t *M);
  /* Sets {M} to the transpose {A^t} of matrix {A}. */

void r2x2_get_row(r2x2_t *A, int32_t i, r2_t *x);
void r2x2_set_row(r2x2_t *A, int32_t i, r2_t *x);
  /* These two procedures copy row {i} of matrix {A} to and from vector {x}, respectively. */

void r2x2_get_col(r2x2_t *A, int32_t j, r2_t *x);
void r2x2_set_col(r2x2_t *A, int32_t j, r2_t *x);
  /* These two procedures copy column {j} of matrix {A} to and from vector {x}, respectively. */

void r2x2_map_row (r2_t *x, r2x2_t *A, r2_t *r);
  /* Sets {r} to the product of row vector {x} by matrix {A}. */

void r2x2_map_col (r2x2_t *A, r2_t *x, r2_t *r);
  /* Sets {r} to the product of matrix {A} by column vector {x}. */
  
void r2x2_add (r2x2_t *A, r2x2_t *B, r2x2_t *M);
  /* Sets {M = A + B}. */

void r2x2_sub (r2x2_t *A, r2x2_t *B, r2x2_t *M);
  /* Sets {M = A - B}. */

void r2x2_neg (r2x2_t *A, r2x2_t *M);
  /* Sets {M} to {-A}. */

void r2x2_scale (double s, r2x2_t *A, r2x2_t *M);
  /* Sets {M} to the product of scalar {s} and matrix {A}. */

void r2x2_mix (double s, r2x2_t *A, double t, r2x2_t *B, r2x2_t *M);
  /* Sets {M := s * A + t * B}. */

void r2x2_mul (r2x2_t *A, r2x2_t *B, r2x2_t *M);
  /* Sets {M} to the product of matrices {A} and {B}. */

void r2x2_mul_tr (r2x2_t *A, r2x2_t *B, r2x2_t *M);
  /* Computes the matrix product {M = A * B^t}. (In other words, sets
    {M[i,j]} to the dot product of row {i} of {A} and row {j} of {B}.) */

void r2x2_tr_mul (r2x2_t *A, r2x2_t *B, r2x2_t *M);
  /* Computes the matrix product {M = A^t * B}. (In other words, sets
    {M[i,j]} to the dot product of column {i} of {A} and column {j} of {B}.) */

double r2x2_det (r2x2_t *A);
  /* Returns the determinant of matrix {A}. */

void r2x2_adj (r2x2_t *A, r2x2_t *M);
  /* Sets {M} to the adjoint of matrix {A}. */

void r2x2_inv (r2x2_t *A, r2x2_t *M);
  /* Sets {M} to the inverse of matrix {A}. */

void r2x2_sqrt(r2x2_t *A, r2x2_t *M);
  /* Sets {M} to the square root of {A}. Bombs if the square root
    is non-real or does not exist. */

double r2x2_norm_sqr(r2x2_t* A);
  /* Squared Frobenius norm of {A}, i.e. sum of squares of elements */

double r2x2_norm(r2x2_t* A);
  /* Frobenius norm of {A}, i.e. square root of sum of squares of elements */

double r2x2_normalize(r2x2_t *A);
  /* Divides the matrix {A} by its {r2x2_norm}. 
    If that norm is zero, the matrix is filled with {NAN}.
    Returns the norm. */

double r2x2_mod_norm_sqr (r2x2_t *A);
  /* Returns the square of the Frobenius norm of {A-I} */

bool_t r2x2_is_unif_scaling(r2x2_t *M, double s);
  /* TRUE iff {M} is a diagonal matrix with all diagonal
    elements equal to {s}. */

void r2x2_from_rows(r2_t *a, r2_t *b, r2x2_t *M);
  /* Sets {M} to the matrix whose rows are the vectors {a,b}. */

void r2x2_from_cols(r2_t *a, r2_t *b, r2x2_t *M);
  /* Sets {M} to the matrix whose columns are the vectors {a,b}. */

void r2x2_print (FILE *f, r2x2_t *A);
  /* Prints matrix {A} to file {f}, with default format. */

void r2x2_gen_print
  ( FILE *f, r2x2_t *A,
    char *fmt, 
    char *olp, char *osep, char *orp, /* Outer delimiters. */
    char *ilp, char *isep, char *irp  /* Inner delimiters. */
  );
  /* Prints matrix {A} to file {f}, using {fmt} for each element.
    The matrix is bounded by {olp} and {orp}, and rows are separated
    by {osep}.  Each row is bounded by {ilp} and {irp}, and elements
    are separated by {isep}. Defaults are provided for any of these
    strings which are NULL. */

/* OPERATIONS SPECIFIC TO 2x2 MATRICES */

void r2x2_rot90(r2x2_t *M);
  /* Stores in {M} the matrix {[[00,+1],[-1,00]]}, so that
    {r2x2_map_row(x,M,x)} rotates {x} by 90 degrees counterclockwise. */

void r2x2_rot_and_scale(r2_t *p, r2x2_t *M);
  /* Stores in {M} the similarity (rotation and uniform scaling) matrix
    that keeps the origin fixed and maps {(1,0)} to {p}. */

void r2x2_moments(r2x2_t *A, r2_t *e);
  /* Stores in {e} the eigenvalues of {A*A^t} where {A^t} is the 
    transpose of {A}.  The eigenvalues are non-negative and
    sorted in decreasing order. They are the max and min 
    values of {|u*A|^2/|u|^2} for any nonzero vector {u}. */

void r2x2_sym_eigen(r2x2_t *A, r2_t *e, r2x2_t *M);
  /* Computes the eigenvalues {e.c[0],e.c[1]) of the symmetric matrix
    {A}. The eigenvalues will be sorted in decreasing order of SIGNED
    value.
    
    Actually, the procedure computes the eigenvalues of the matrix
    {(A+A')/2}, where {A'} is the transpose of {A}.  That's always
    symmetric, and equal to {A} if {A} is symmetric.
    
    If {M} is not NULL, the procedure also stores into the rows of {M}
    the corresponding unit-length eigenvectors. In that case, {M} will
    be a rotation matrix (i.e. {M'*M = M*M'= IDENT()}, {det(M)=+1})
    such that {A = M'*DIAG(e)*M} (or {M*A*M' = DIAG(e)}), where
    {DIAG(e)} is the diagonal matrix with the elements of {e} along
    the diagonal. */

void r2x2_from_point_pairs(r2_vec_t *p1, r2_t *bar1, r2_vec_t *p2, r2_t *bar2, r2x2_t *M);
  /* Returns in {*M} the linear map that takes the points {p1} as close as
    possible to the points {p2}, in the sense of minimizing the mean
    squared error.
    
    If {bar1} is not NULL, the point {*bar1} is subtracted from every
    point in {p1}.  Ditto for {bar2} and {p2}.
    
    The lists {p1 and {p2} must have the same length, and {p1} (after
    subtracting {bar1}, if not null) must contain at least two
    linearly independent points; otherwise the result may
    contain infinite or {NAN} elements. */

#endif

