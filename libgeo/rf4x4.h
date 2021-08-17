/* rf4x4.h --- 4x4 matrices and operations on them (single-precision version) */
/* Last edited on 2021-08-17 08:49:56 by stolfi */

#ifndef rf4x4_H
#define rf4x4_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <rf4.h>

typedef struct { float c[4][4]; } rf4x4_t;

/* STANDARD OPERATIONS */

rf4x4_t rf4x4_zero(void);
  /* Returns the null matrix. */

rf4x4_t rf4x4_ident(void);
  /* Returns the identity matrix. */

rf4x4_t rf4x4_transp (rf4x4_t *A);
  /* Returns the transpose {A^t} of matrix {A} */

rf4_t rf4x4_get_row(rf4x4_t *A, int32_t i);
  /* Extracts row {i} of matrix {A} as a vector. */

void rf4x4_set_row(rf4x4_t *A, int32_t i, rf4_t *x);
  /* Sets row {i} of matrix {A} to the vector {x}. */

rf4_t rf4x4_get_col(rf4x4_t *A, int32_t j);
  /* Extracts column {j} of matrix {A} as a vector. */

void rf4x4_set_col(rf4x4_t *A, int32_t j, rf4_t *x);
  /* Sets column {j} of matrix {A} to vector {x}. */

rf4_t rf4x4_map_row (rf4_t *x, rf4x4_t *A);
  /* Returns the product of row vector {x} by matrix {A} */

rf4_t rf4x4_map_col (rf4x4_t *A, rf4_t *x);
  /* Retruns the product of matrix {A} by column vector {x} */

rf4x4_t rf4x4_scale (double s, rf4x4_t *A);
  /* Returns the product of scalar {s} and matrix {A}. */
 
rf4x4_t rf4x4_mul (rf4x4_t *A, rf4x4_t *B);
  /* Returns the product of matrices {A} and {B} */

rf4x4_t rf4x4_mul_tr (rf4x4_t *A, rf4x4_t *B);
  /* Returns the matrix product {M = A * B^t}. (In other words, sets
    {M[i,j]} to the dot product of row {i} of {A} and row {j} of {B}.) */

double rf4x4_det (rf4x4_t *A);
  /* Returns the determinant of matrix {A} */

rf4x4_t rf4x4_adj (rf4x4_t *A);
  /* Returns the adjoint of matrix {A} */

rf4x4_t rf4x4_inv (rf4x4_t *A);
  /* Returns the inverse of matrix {A} */

double rf4x4_norm_sqr(rf4x4_t* A);
  /* Squared Frobenius norm of {A}, i.e. sum of squares of elements */

double rf4x4_norm(rf4x4_t* A);
  /* Frobenius norm of {A}, i.e. square root of sum of squares of elements */

double rf4x4_mod_norm_sqr (rf4x4_t *A);
  /* Returns the square of the Frobenius norm of {A-I} */

bool_t rf4x4_is_unif_scaling(rf4x4_t *M, double s);
  /* TRUE iff {M} is a diagonal matrix with all diagonal
    elements equal to {s}. */

rf4x4_t rf4x4_from_rows(rf4_t *a, rf4_t *b, rf4_t *c, rf4_t *d);
  /* Sets {M} to the matrix whose rows are the vectors {a,b,c,d}. */

rf4x4_t rf4x4_from_cols(rf4_t *a, rf4_t *b, rf4_t *c, rf4_t *d);
  /* Sets {M} to the matrix whose columns are the vectors {a,b,c,d}. */

void rf4x4_print (FILE *f, rf4x4_t *A);
  /* Prints matrix {A} to file {f}, with default format. */

void rf4x4_gen_print
  ( FILE *f, rf4x4_t *A,
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

