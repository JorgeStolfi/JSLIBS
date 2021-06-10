/* r4x4.h --- 4x4 matrices and operations on them */
/* Last edited on 2021-06-09 19:44:44 by jstolfi */

#ifndef r4x4_H
#define r4x4_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <r4.h>

typedef struct { double c[4][4]; } r4x4_t;

/* STANDARD OPERATIONS */

void r4x4_zero(r4x4_t *M);
  /* Stores in {M} the null matrix. */

void r4x4_ident(r4x4_t *M);
  /* Stores in {M} the identity matrix. */

void r4x4_transp (r4x4_t *A, r4x4_t *M);
  /* Sets {M} to the transpose {A^t} of matrix {A} */

void r4x4_get_row(r4x4_t *A, int32_t i, r4_t *x);
void r4x4_set_row(r4x4_t *A, int32_t i, r4_t *x);
  /* These two procedures copy row {i} of matrix {A} to and from vector {x}, respectively. */

void r4x4_get_col(r4x4_t *A, int32_t j, r4_t *x);
void r4x4_set_col(r4x4_t *A, int32_t j, r4_t *x);
  /* These two procedures copy column {j} of matrix {A} to and from vector {x}, respectively. */

void r4x4_map_row (r4_t *x, r4x4_t *A, r4_t *r);
  /* Sets {r} to the product of row vector {x} by matrix {A} */

void r4x4_map_col (r4x4_t *A, r4_t *x, r4_t *r);
  /* Sets {r} to the product of matrix {A} by column vector {x} */

void r4x4_scale (double s, r4x4_t *A, r4x4_t *M);
  /* Sets {M} to the product of scalar {s} and matrix {A}. */
 
void r4x4_mul (r4x4_t *A, r4x4_t *B, r4x4_t *M);
  /* Sets {M} to the product of matrices {A} and {B} */

void r4x4_mul_tr (r4x4_t *A, r4x4_t *B, r4x4_t *M);
  /* Computes the matrix product {M = A * B^t}. (In other words, sets
    {M[i,j]} to the dot product of row {i} of {A} and row {j} of {B}.) */

double r4x4_det (r4x4_t *A);
  /* Returns the determinant of matrix {A} */

void r4x4_adj (r4x4_t *A, r4x4_t *M);
  /* Sets {M} to the adjoint of matrix {A} */

void r4x4_inv (r4x4_t *A, r4x4_t *M);
  /* Sets {M} to the inverse of matrix {A} */

double r4x4_norm_sqr(r4x4_t* A);
  /* Squared Frobenius norm of {A}, i.e. sum of squares of elements */

double r4x4_norm(r4x4_t* A);
  /* Frobenius norm of {A}, i.e. square root of sum of squares of elements */

double r4x4_mod_norm_sqr (r4x4_t *A);
  /* Returns the square of the Frobenius norm of {A-I} */

bool_t r4x4_is_unif_scaling(r4x4_t *M, double s);
  /* TRUE iff {M} is a diagonal matrix with all diagonal
    elements equal to {s}. */

void r4x4_from_rows(r4_t *a, r4_t *b, r4_t *c, r4_t *d, r4x4_t *M);
  /* Sets {M} to the matrix whose rows are the vectors {a,b,c,d}. */

void r4x4_from_cols(r4_t *a, r4_t *b, r4_t *c, r4_t *d, r4x4_t *M);
  /* Sets {M} to the matrix whose columns are the vectors {a,b,c,d}. */

void r4x4_print (FILE *f, r4x4_t *A);
  /* Prints matrix {A} to file {f}, with default format. */

void r4x4_gen_print
  ( FILE *f, r4x4_t *A,
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

