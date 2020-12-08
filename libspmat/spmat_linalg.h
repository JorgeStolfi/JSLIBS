#ifndef spmat_linalg_H
#define spmat_linalg_H
/* Linear algebra operations on generic sparse matrices. */

#define spmat_linalg_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2009-08-31 21:42:26 by stolfi */

/* LINEAR ALGEBRA OPERATIONS

  This interface provides some basic linear algebra functions
  on generic sparse matrices.  For example, if {dspmat_t}
  was defined as a sparse matrix with {double} elements,
  one could write:

    ------------------------------------------------------------
    //Create an {n × n} identity matrix {I}:
    dspmat_t I = dspmat_new(n,n,n);
    dspmat_identity(&I,n);

    //Set {R = M - 2*I}:
    dspmat_t R = dspmat_new(n,n,0);
    dspmat_mix(1.0, &M, -2.0, &I, &R);

    //Set {S = R*R}:
    dspmat_t S = dspmat_new(n,n,0);
    dspmat_mul(&R, &R, &S);

    //Write the result to {stdout}:
    dspmat_write(stdout, &S);
    ------------------------------------------------------------ */
  /* To obtain these linear algebra functions, one should use

    ------------------------------------------------------------
    / * Declarations * /
    #include <spmat_linalg.h>

    spmat_typedef(dspmat_t, dspmat, double);

    spmat_linalg_def(dspmat_t, dspmat, double);
    ------------------------------------------------------------ */
  /*   ------------------------------------------------------------
    / * Implementations  * /
    #include <spmat_linalg.h>

    spmat_impl(dspmat_t, dspmat, double);

    #define dspmat_trivial_elem() (0.0)
    #define dspmat_elem_is_trivial(X) ((X)==0.0)

    #define dspmat_elem_zero (0.0)
    #define dspmat_elem_one (1.0)
    #define dspmat_elem_add(X,Y) ((X)+(Y))
    #define dspmat_elem_mul(X,Y) ((X)*(Y))

    spmat_impl(dspmat_t, dspmat, double);
    spmat_linalg_impl(dspmat_t, dspmat, double);
    ------------------------------------------------------------

 */

#include <stdlib.h>
#include <stdint.h>

#include <spmat.h>

#define spmat_linalg_def(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_LINALG_DEF(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/*
  This macro expands into prototype declarations
  of the basic linear algebra procedures:

    ------------------------------------------------------------
    PREFIX##_identity
    PREFIX##_mix
    PREFIX##_mul
    PREFIX##_map_col
    PREFIX##_map_row
    ------------------------------------------------------------
  
  They are appropriate when the addition and multiplication of
  {ELEM_TYPE} values are defined. */

#define spmat_linalg_impl(MATRIX_TYPE,PREFIX,ELEM_TYPE) \
  spmat_LINALG_IMPL(MATRIX_TYPE,PREFIX,ELEM_TYPE)
/*
  This macro expands into the implementations of the linear algebra
  procedures declared by {spmat_linalg_def}.

  BEFORE calling this macro, the client must declare the following
  macros or procedures:

    ------------------------------------------------------------
    #define PREFIX##_elem_add(X,Y) (...)
       / * The addition of two {ELEM_TYPE} values {X,Y}. * /

    #define PREFIX##_elem_mul(X,Y) (...)
       / * The multiplication of two {ELEM_TYPE} values {X,Y}. * /

    #define PREFIX##_elem_one (...)
      / * The {ELEM_TYPE} value that is the multiplicative unit. * /
    ------------------------------------------------------------
  
  This interface assumes that {PREFIX##_trivial_elem} (the element
  value which is not stored explicitly in the matrix) is the
  zero value for {PREFIX##_elem_add}. */

/* ======================================================================

  The remainder of this interface describes the functions provided by
  the macros {spmat_linalg_def} and {spmat_linalg_impl}. */

/* IDENTITY MATRIX*/

/* ------------------------------------------------------------
void PREFIX##_identity(MATRIX_TYPE *M, spmat_size_t size);
------------------------------------------------------------ */
  /* Sets the matrix {M} to the identity matrix with
    {size} rows and {size} columns.*/

/* ------------------------------------------------------------
void PREFIX##_mix
  ( ELEM_TYPE a,
    MATRIX_TYPE *A,
    ELEM_TYPE b,
    MATRIX_TYPE *B,
    MATRIX_TYPE *C
  );
------------------------------------------------------------ */
  /* Stores into {C} the linear combination {a*A + b*B}.

    The entries of {A} MUST BE SORTED by increasing row index, then by
    increasing column index within each row (as produced by
    {PREFIX##_sort_entries(&A, +2, +1)}); and ditto for {B}. The entries
    of {C} will be in the same order.  If either matrix contains
    two or more entries with the same indices, the result
    may have repeated entries, too.
    
    The matrices {A} and {B} may be the same, but the matrix {C} must
    be storage-disjoint from both.

    The matrices {A} and {B} must also have the same dimensions
    {A.rows==B.rows} and {A.cols==B.cols}, which will be assigned to
    {C.rows} and {C.cols}. The vector {C.e} will be re-alocated, if
    necessary, to hold exactly the non-trivial entries of the result.*/

/* ------------------------------------------------------------
void PREFIX##_mul(MATRIX_TYPE *A, MATRIX_TYPE *B, MATRIX_TYPE *C);
------------------------------------------------------------ */
  /* Stores into {C} the matrix product {A*B}.

    The entries of {A} (only) MUST BE SORTED by increasing row index
    (as obtained with {PREFIX##_sort_entries(&A, +1, 0)}.
    The entries of {C} will be sorted by increasing row index,
    then by increasing column index within each row.  If either matrix contains
    two or more entries with the same indices, the result
    may have repeated entries, too.
    
    The matrices {A} and {B} may be the same, but the matrix {C} must
    be storage-dsjoint from both.

    The procedure requires {A.cols == B.rows}, and will set
    {R.rows=A.rows} and {R.cols=B.cols}. The vector {R.e} will be
    re-alocated if necessary to hold exactly the non-trivial entries
    of the result. */

/* MATRIX-VECTOR MULTIPLICATION*/

/* ------------------------------------------------------------
void PREFIX##_map_col
  ( MATRIX_TYPE *M,
    ELEM_TYPE a[],
    spmat_size_t na,
    ELEM_TYPE b[],
    spmat_size_t nb
  );
------------------------------------------------------------ */
  /* This macro declares the function {PREFIX##_map_col}. The call
    {PREFIX##_map_col(&M, a, na, b, nb)} will compute
    the product of the matrix {M} by the column vector
    {a[0..na-1]} and store the result into {b[0..nb-1]}. 
    
    If {M} contains two or more entries with the same indices, the
    corresponding contributions to {b} are added together.

    The procedure requires {na == M.cols} and {nb == M.rows}. The
    vectors {a} and {b} must be storage-disjoint. The entries of {M}
    need not be sorted. */

/* ------------------------------------------------------------
void PREFIX##_map_row
  ( ELEM_TYPE a[],
    spmat_size_t na,
    MATRIX_TYPE *M,
    ELEM_TYPE b[],
    spmat_size_t nb
  );
------------------------------------------------------------ */
  /* This macro declares the function {PREFIX##_map_row}. The call
    {PREFIX##_map_row(a, na, &M, b, nb)} will compute
    the product of the row vector
    {a[0..na-1]} by the matrix {M}, and store the result into {b[0..nb-1]}.
    
    If {M} contains two or more entries with the same indices, the
    corresponding contributions to {b} are added together.

    The procedure requires {na == M.rows} and {nb == M.cols}. The
    vectors {a} and {b} must be storage-disjoint. The entries of {M}
    need not be sorted. */

#include <spmat_linalg_def.h>

#endif
