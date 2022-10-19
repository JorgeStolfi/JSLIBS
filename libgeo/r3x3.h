/* r3x3.h --- 3x3 matrices and operations on them */
/* Last edited on 2022-01-05 14:15:54 by stolfi */

#ifndef r3x3_H
#define r3x3_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <r3.h>

typedef struct { double c[3][3]; } r3x3_t;

/* STANDARD OPERATIONS */

void r3x3_zero(r3x3_t *M);
  /* Stores in {M} the null matrix. */

void r3x3_ident(r3x3_t *M);
  /* Stores in {M} the identity matrix. */

void r3x3_transp (r3x3_t *A, r3x3_t *M);
  /* Sets {M} to the transpose {A^t} of matrix {A} */

void r3x3_get_row(r3x3_t *A, int32_t i, r3_t *x);
void r3x3_set_row(r3x3_t *A, int32_t i, r3_t *x);
  /* These two procedures copy row {i} of matrix {A} to and from vector {x}, respectively. */

void r3x3_get_col(r3x3_t *A, int32_t j, r3_t *x);
void r3x3_set_col(r3x3_t *A, int32_t j, r3_t *x);
  /* These two procedures copy column {j} of matrix {A} to and from vector {x}, respectively. */

void r3x3_map_row (r3_t *x, r3x3_t *A, r3_t *r);
  /* Sets {r} to the product of row vector {x} by matrix {A} */

void r3x3_map_col (r3x3_t *A, r3_t *x, r3_t *r);
  /* Sets {r} to the product of matrix {A} by column vector {x} */
  
void r3x3_add (r3x3_t *A, r3x3_t *B, r3x3_t *M);
  /* Sets {M = A + B}. */

void r3x3_sub (r3x3_t *A, r3x3_t *B, r3x3_t *M);
  /* Sets {M = A - B}. */

void r3x3_neg (r3x3_t *A, r3x3_t *M);
  /* Sets {M} to {-A}. */

void r3x3_scale (double s, r3x3_t *A, r3x3_t *M);
  /* Sets {M} to the product of scalar {s} and matrix {A}. */

void r3x3_mix (double s, r3x3_t *A, double t, r3x3_t *B, r3x3_t *M);
  /* Sets {M} to the linear combination {s*A + t*B}. */

void r3x3_mul (r3x3_t *A, r3x3_t *B, r3x3_t *M);
  /* Sets {M} to the product of matrices {A} and {B} */

void r3x3_mul_tr (r3x3_t *A, r3x3_t *B, r3x3_t *M);
  /* Computes the matrix product {M = A * B^t}. (In other words, sets
    {M[i,j]} to the dot product of row {i} of {A} and row {j} of {B}.) */

double r3x3_det (r3x3_t *A);
  /* Returns the determinant of matrix {A} */

void r3x3_adj (r3x3_t *A, r3x3_t *M);
  /* Sets {M} to the adjoint of matrix {A} */

void r3x3_inv (r3x3_t *A, r3x3_t *M);
  /* Sets {M} to the inverse of matrix {A} */

double r3x3_norm_sqr(r3x3_t* A);
  /* Squared Frobenius norm of {A}, i.e. sum of squares of elements */

double r3x3_norm(r3x3_t* A);
  /* Frobenius norm of {A}, i.e. square root of sum of squares of elements */

double r3x3_mod_norm_sqr (r3x3_t *A);
  /* Returns the square of the Frobenius norm of {A-I} */

void r3x3_diff_sqr(r3x3_t *A, r3x3_t *B, r3x3_t *R, double *dabs2P, double *drel2P);
  /* Returns in {*dabs2P} and {*drel2P} (if not NULL)
    the total squared discrepancy between the matrices {*A} and 
    {*B}, respectively absolute and relative to the matrix {*R}.  
    Namely 
    
      { dabs2 = SUM { (Ae[s] - Be[s])^2 } }
    
      { drel2 = SUM { ((Ae[s] - Be[s])/Re[s])^2 } }
    
    where {Ae[s],Be[s],Re[s]} are all corresponding elements of
    {*A,*B,*C}. excluding those where {Re[s]} is zero. */

bool_t r3x3_is_unif_scaling(r3x3_t *A, double s);
  /* TRUE iff {A} is a diagonal matrix with all diagonal
    elements equal to {s}. */

void r3x3_from_rows(r3_t *a, r3_t *b, r3_t *c, r3x3_t *M);
  /* Sets {M} to the matrix whose rows are the vectors {a,b,c}. */

void r3x3_from_cols(r3_t *a, r3_t *b, r3_t *c, r3x3_t *M);
  /* Sets {M} to the matrix whose columns are the vectors {a,b,c}. */
  
/* OPERATIONS SPECIFIC TO 3x3 MATRICES */

void r3x3_u_v_rotation(r3_t *u, r3_t *v, r3x3_t *M);
  /* Sets {M} to a rotation matrix, around some axis through the
    origin, that takes the unit vector {u} to to the unit vector {v}
    by the shortest route.
    
    If {u} is equal to {v}, then {M} will be the identity matrix.
    If {u} is opposite to {v}, {M} will be a 180 degree rotation around a
    random axis through the origin that is orthogonal to both. */

/* I/O */

void r3x3_print (FILE *f, r3x3_t *A);
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

#endif

