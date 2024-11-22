/* sym_eigen.h -- eigenvalues and eigenvectors of a symmetric matrix */
/* Last edited on 2024-11-20 15:03:30 by stolfi */

#ifndef sym_eigen_H
#define sym_eigen_H

#include <stdint.h>

void sym_eigen_tridiagonalize
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

       {A[0..n*n-1]} is the input symmetric matrix, stored by rows.  
         (The procedure will only use the lower triangular
         part of {A}, including the diagonal.)

    On output,

       {d[0..n-1]} are the diagonal elements of the tridiagonal
         matrix {T}.  (I.e. {d[i] == T[i,i]}, for {i=0..n-1}.) 

       {e[1..n-1]} are the subdiagonal elements of {T}.
         (I.e. {e[i] == T[i,i-1]=T[i-1,i]}, for {i=1..n-1}.) 
         Element {e[0]} is set to zero.

       {R[0..n*n-1} contains the
         orthogonal transformation matrix used for the reduction, i.e.
         the matrix such that {R*A*(R^t) == T}.

    The matrix {R} may be the same as {A}, which is then overwritten
    with the trabsformation matrix. Otherwise, {A} is unaltered.
  */

void sym_eigen_trid_eigen
  ( uint32_t n,
    double *d,
    double *e,
    double *R,
    uint32_t *p,
    uint32_t absrt
  );
  /* Finds the right eigenvalues and eigenvectors of a symmetric
    tridiagonal matrix {T} by the QL method. 
    
    If the matrix {T} vas derived from a full symmetric matrix {A} by
    an orthogonal similarity map {R} (see {sym_eigen_tridiagonalize}), this
    procedure will find the eigenvectors of {A}.

    On input,

       {n} is the size (rows and columns) of the matrices {T} and {R}.

       {d} must contain the diagonal elements of the input matrix {T}.

       {e} must contain the subdiagonal elements of the input matrix {T}
         in its last {n-1} positions (i.e. {e[i] == T[i,i-1]},
         for {i=1..n-1}).  The value of {e[0]} is ignored.

       {R} must contain the orthogonal similarity matrix which
         transformed the original matrix {A} into {T}, i.e. the matrix
         such that {R*A*(R^t) == T}. (To obtain the eigenvectors of
         {T} itself, {R} must be set to the identity matrix.)
         
       {absrt} indicates the desired ordering of the eigenvalues: 
         1 means by absolute value, 0 by signed value. 

     On output,

       {*p} is the number of eigenvalues actually computed.
         (If less than {n}, then the computation of some
         some eigenvalue failed to converge after 30 iterations).

       {d[0..*p-1]} contains the eigenvalues that could be 
         determined, in ascending order of absolute value
         (if {absrt == 1}) or signed value (if {absrt == 0}).

       {e} has been destroyed.

       {R[0..*p-1,0..n-1]} are orthonormal eigenvectors of the matrix
         {A}, corresponding to the eigenvalues {d[0..*p-1]}, stored as
         rows. If {S} is {R} truncated to its first {p} rows, then
         {S*A*(S^t) == diag(d[0..p-1])}. */

#endif
