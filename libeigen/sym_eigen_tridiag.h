/* sym_eigen_tridiag.h -- tridiagonalization of a symmetric matrix */
/* Last edited on 2024-12-05 11:10:58 by stolfi */

#ifndef sym_eigen_tridiag_H
#define sym_eigen_tridiag_H

#include <stdint.h>

void sym_eigen_tridiag
  ( uint32_t n,
    double A[],
    double d[],
    double e[],
    double R[]
  );
  /* Reduces the real symmetric matrix {A} to a symmetric tridiagonal
    matrix {T}, using an orthogonal similarity transformation {R}.

    On input,

       {n} is the size (rows and columns) of the matrices {A} and {R}.
         Must be positive.

       {A[0..n*n-1]} is the input symmetric matrix, stored by rows.  
         (The procedure will only use the lower triangular
         part of {A}, including the diagonal.)

    On output,

       {d[0..n-1]} are the diagonal elements of the tridiagonal
         matrix {T}.  (I.e. {d[i] == T[i,i]}, for {i=0..n-1}.) 

       {e[1..n-1]} are the subdiagonal elements of {T}.
         (I.e. {e[i] == T[i,i-1]=T[i-1,i]}, for {i=1..n-1}.) 
         Element {e[0]} is set to zero.

       {R[0..n*n-1} contains the orthogonal transformation 
         matrix used for the reduction, i.e. the matrix such
         that {R*A*(R^t) == T}.

    The matrix {R} may be {NULL}, in which case the relation between {A}
    and {T} is lost.  Also {R} may be the same as {A}, which is then overwritten
    with the trabsformation matrix. Otherwise, {A} is unaltered. */

#endif
