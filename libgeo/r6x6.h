/* r6x6.h --- 4x4 matrices and operations on them */
/* Last edited on 2024-11-20 13:01:29 by stolfi */

#ifndef r6x6_H
#define r6x6_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <sign.h>
#include <r6.h>

typedef struct r6x6_t { double c[6][6]; } r6x6_t;

/* STANDARD OPERATIONS */

void r6x6_zero(r6x6_t *M);
  /* Stores in {M} the null matrix. */

void r6x6_ident(r6x6_t *M);
  /* Stores in {M} the identity matrix. */

void r6x6_throw(r6x6_t *M, sign_t sgn);
  /* Fills {M} with random elements in the range {[-1 _ +1]}.
    The {sgn} must be {-1}, 0, or {+1}. If it is not zero,
    the matrix will have a nonzero determinant of that sign.
    If {sgn} is zero, the determinant may have any sign,
    including zero. */

void r6x6_transp (r6x6_t *A, r6x6_t *M);
  /* Sets {M} to the transpose {A^t} of matrix {A} */

void r6x6_get_row(r6x6_t *A, uint32_t i, r6_t *x);
void r6x6_set_row(r6x6_t *A, uint32_t i, r6_t *x);
  /* These two procedures copy row {i} of matrix {A} to and from vector {x}, respectively. */

void r6x6_get_col(r6x6_t *A, uint32_t j, r6_t *x);
void r6x6_set_col(r6x6_t *A, uint32_t j, r6_t *x);
  /* These two procedures copy column {j} of matrix {A} to and from vector {x}, respectively. */

void r6x6_map_row (r6_t *x, r6x6_t *A, r6_t *r);
  /* Sets {r} to the product of row vector {x} by matrix {A} */

void r6x6_map_col (r6x6_t *A, r6_t *x, r6_t *r);
  /* Sets {r} to the product of matrix {A} by column vector {x} */

void r6x6_scale (double s, r6x6_t *A, r6x6_t *M);
  /* Sets {M} to the product of scalar {s} and matrix {A}. */

void r6x6_mul (r6x6_t *A, r6x6_t *B, r6x6_t *M);
  /* Sets {M} to the product of matrices {A} and {B} */

void r6x6_mul_tr (r6x6_t *A, r6x6_t *B, r6x6_t *M);
  /* Computes the matrix product {M = A * B^t}. (In other words, sets
    {M[i,j]} to the dot product of row {i} of {A} and row {j} of {B}.) */

double r6x6_det (r6x6_t *A);
  /* Returns the determinant of matrix {A} */

/* !!! Uncomment when we get {rmxn_adj}: !!!
void r6x6_adj (r6x6_t *A, r6x6_t *M);
*/
  /* Sets {M} to the adjoint of matrix {A} */

void r6x6_inv (r6x6_t *A, r6x6_t *M);
  /* Sets {M} to the inverse of matrix {A} */

double r6x6_norm_sqr(r6x6_t* A);
  /* Squared Frobenius norm of {A}, i.e. sum of squares of elements */

double r6x6_norm(r6x6_t* A);
  /* Frobenius norm of {A}, i.e. square root of sum of squares of elements */

double r6x6_normalize(r6x6_t *A);
  /* Divides the matrix {A} by its {r6x6_norm}. 
    If that norm is zero, the matrix is filled with {NAN}.
    Returns the norm. */
 
double r6x6_mod_norm_sqr (r6x6_t *A);
  /* Returns the square of the Frobenius norm of {A-I} */
  
bool_t r6x6_is_unif_scaling(r6x6_t *M, double s);
  /* TRUE iff {M} is a diagonal matrix with all diagonal
    elements equal to {s}. */

void r6x6_from_rows(r6_t *a, r6_t *b, r6_t *c, r6_t *d, r6_t *e, r6_t *f, r6x6_t *M);
  /* Sets {M} to the matrix whose rows are the vectors {a,b,c,d,e,f}. */

void r6x6_from_cols(r6_t *a, r6_t *b, r6_t *c, r6_t *d, r6_t *e, r6_t *f, r6x6_t *M);
  /* Sets {M} to the matrix whose columns are the vectors {a,b,c,d,e,f}. */
  
void r6x6_print (FILE *f, r6x6_t *A);
  /* Prints matrix {A} to file {f}, with default format. */

void r6x6_gen_print
  ( FILE *f, r6x6_t *A,
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

