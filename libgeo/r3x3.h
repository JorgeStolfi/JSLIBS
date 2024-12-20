/* r3x3.h --- 3x3 matrices and operations on them */
/* Last edited on 2024-12-05 10:28:04 by stolfi */

#ifndef r3x3_H
#define r3x3_H

#include <stdio.h>
#include <stdint.h>

#include <sign.h>
#include <r3.h>

typedef struct r3x3_t { double c[3][3]; } r3x3_t;

/* STANDARD OPERATIONS */

void r3x3_zero(r3x3_t *M);
  /* Stores in {M} the null matrix. */

void r3x3_ident(r3x3_t *M);
  /* Stores in {M} the identity matrix. */

void r3x3_throw(r3x3_t *M, sign_t sgn);
  /* Fills {M} with random elements in the range {[-1 _ +1]}.
    The {sgn} must be {-1}, 0, or {+1}. If it is not zero,
    the matrix will have a nonzero determinant of that sign.
    If {sgn} is zero, the determinant may have any sign,
    including zero. */

void r3x3_throw_rotation(r3x3_t *M);
  /* Fills {M} with a random orthonormal matrix with 
    positive determinant ({+1} apart from roundoff errors).

    Either pre- or post- multiplication of a 3-vector by the matrix will
    perform a rotation of {\RR^3} by a random angle about some random
    axis through the origin. The identity matrix is a (highly unlikely)
    special case, where the angle is zero and the axis is
    indeterminate. */

void r3x3_transp(r3x3_t *A, r3x3_t *M);
  /* Sets {M} to the transpose {A^t} of matrix {A}. */

void r3x3_get_row(r3x3_t *A, uint32_t i, r3_t *x);
void r3x3_set_row(r3x3_t *A, uint32_t i, r3_t *x);
  /* These two procedures copy row {i} of matrix {A} to and from vector {x}, respectively. */

void r3x3_get_col(r3x3_t *A, uint32_t j, r3_t *x);
void r3x3_set_col(r3x3_t *A, uint32_t j, r3_t *x);
  /* These two procedures copy column {j} of matrix {A} to and from vector {x}, respectively. */

void r3x3_map_row(r3_t *x, r3x3_t *A, r3_t *r);
  /* Sets {r} to the product of row vector {x} by matrix {A}. */

void r3x3_map_col(r3x3_t *A, r3_t *x, r3_t *r);
  /* Sets {r} to the product of matrix {A} by column vector {x}. */
  
void r3x3_add(r3x3_t *A, r3x3_t *B, r3x3_t *M);
  /* Sets {M} to the matrix sum {A + B}. */

void r3x3_sub(r3x3_t *A, r3x3_t *B, r3x3_t *M);
  /* Sets {M} to the matrix difference {A - B}. */

void r3x3_neg(r3x3_t *A, r3x3_t *M);
  /* Sets {M} to the negated matrix {-A}. */

void r3x3_scale(double s, r3x3_t *A, r3x3_t *M);
  /* Sets {M} to the product of scalar {s} and matrix {A}. */

void r3x3_mix(double s, r3x3_t *A, double t, r3x3_t *B, r3x3_t *M);
  /* Sets {M} to the linear combination {s*A + t*B}. */

void r3x3_mul(r3x3_t *A, r3x3_t *B, r3x3_t *M);
  /* Sets {M} to the product of matrices {A} and {B}. */

void r3x3_mul_tr(r3x3_t *A, r3x3_t *B, r3x3_t *M);
  /* Computes the matrix product {M = A * B^t}. (In other words, sets
    {M[i,j]} to the dot product of row {i} of {A} and row {j} of {B}.) */

void r3x3_tr_mul(r3x3_t *A, r3x3_t *B, r3x3_t *M);
  /* Computes the matrix product {M = A^t * B}. (In other words, sets
    {M[i,j]} to the dot product of column {i} of {A} and column {j} of {B}.) */

double r3x3_det(r3x3_t *A);
  /* Returns the determinant of matrix {A}. */

void r3x3_adj(r3x3_t *A, r3x3_t *M);
  /* Sets {M} to the adjoint of matrix {A}. */

void r3x3_inv(r3x3_t *A, r3x3_t *M);
  /* Sets {M} to the inverse of matrix {A}. */

double r3x3_norm_sqr(r3x3_t *A);
  /* Squared Frobenius norm of {A}, i.e. sum of squares of elements. */

double r3x3_norm(r3x3_t *A);
  /* Frobenius norm of {A}, i.e. square root of sum of squares of elements. */

double r3x3_normalize(r3x3_t *A);
  /* Divides the matrix {A} by its {r3x3_norm}. 
    If that norm is zero, the matrix is filled with {NAN}.
    Returns the norm. */

double r3x3_mod_norm_sqr(r3x3_t *A);
  /* Returns the square of the Frobenius norm of {A-I}. */

double r3x3_L_inf_norm(r3x3_t *A);
  /* Infinity norm of {A}, i.e. max absolute element value. */

double r3x3_L_inf_normalize(r3x3_t *A);
  /* Divides the matrix {A} by its {r3x3_L_inf_norm}. 
    If that norm is zero, the matrix is filled with {NAN}.
    Returns the norm. */

void r3x3_diff_sqr(r3x3_t *A, r3x3_t *B, r3x3_t *R, double *dabs2P, double *drel2P);
  /* Returns in {*dabs2P} and {*drel2P} (if not NULL)
    the total squared discrepancy between the matrices {*A} and 
    {*B}, respectively absolute and relative to the matrix {*R}.  
    Namely 
    
      { dabs2 = SUM { (Ae[s] - Be[s])^2 } }
    
      { drel2 = SUM { ((Ae[s] - Be[s])/Re[s])^2 } }
    
    where {Ae[s],Be[s],Re[s]} are all corresponding elements of
    {*A,*B,*C}. excluding those where {Re[s]} is zero. */

bool_t r3x3_is_unif_scaling(r3x3_t *M, double s, double tol);
  /* TRUE iff {M} is a diagonal matrix with all diagonal
    elements equal to {s}, apart from absolute errors of at most tol. */

void r3x3_from_rows(r3_t *a, r3_t *b, r3_t *c, r3x3_t *M);
  /* Sets {M} to the matrix whose rows are the vectors {a,b,c}. */

void r3x3_from_cols(r3_t *a, r3_t *b, r3_t *c, r3x3_t *M);
  /* Sets {M} to the matrix whose columns are the vectors {a,b,c}. */

/* I/O */

void r3x3_print(FILE *f, r3x3_t *A);
  /* Prints matrix {A} to file {f}, with default format. */

void r3x3_gen_print
  ( FILE *f, r3x3_t *A,
    char *fmt, 
    char *olp, char *osep, char *orp, /* Outer delimiters. */
    char *ilp, char *isep, char *irp  /* Inner delimiters. */
  );
  /* Prints matrix {A} to file {f}, using {fmt} for each element.
    The matrix is bounded by {olp} and {orp}, and rows are separated
    by {osep}.  Each row is bounded by {ilp} and {irp}, and elements
    are separated by {isep}. Defaults are provided for any of these
    strings which are NULL. */
  
/* OPERATIONS SPECIFIC TO 3x3 MATRICES */

void r3x3_u_v_rotation(r3_t *u, r3_t *v, r3x3_t *M);
  /* Sets {M} to a rotation matrix, around some axis through the
    origin, that takes the unit vector {u} to to the unit vector {v}
    by the shortest route.
    
    If {u} is equal to {v}, then {M} will be the identity matrix.
    If {u} is opposite to {v}, {M} will be a 180 degree rotation around a
    random axis through the origin that is orthogonal to both. */

#endif

