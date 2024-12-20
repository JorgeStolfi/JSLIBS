/* sym_eigen.h -- eigenvalues and eigenvectors of a symmetric matrix */
/* Last edited on 2024-12-05 11:10:28 by stolfi */

#ifndef sym_eigen_H
#define sym_eigen_H

#include <stdint.h>

void sym_eigen
  ( uint32_t n,
    double A[],
    double d[],
    double R[],
    uint32_t *nev_P
  );
  /* Finds the right eigenvalues and eigenvectors of a symmetric
    matrix {A} by tridiagonalization followed by the QL method. 

    On input,

       {n} should be the size (rows and columns) of the matrices {A} and {R}.
       
       {A[0..n*n-1]} should be a symmetric {n × n} matrix, stored by rows;
         that is, A{í,j]} should actually be {A[i*n + j]}.

     On output,

       {*nev_P} is the number {nev} of eigenvalues actually computed.
         If less than {n}, then the computation of some
         some eigenvalue failed to converge after 30 iterations.

       {d[0..*p-1]} contains the eigenvalues that could be 
         determined, in ascending order of signed value.

       {R[0..nev-1,0..n-1]} are orthonormal eigenvectors of the matrix
         {A}, corresponding to the eigenvalues {d[0..nev-1]}, stored as
         rows; that is, element {R[i,j]} (coordinate {j} of eigenvector
         {i}) will actually be {R[i*n + j]}. If {S} is {R} truncated to
         its first {nev} rows, then {S*A*(S^t)} will be
         {diag(d[0..nev-1])}.
         
      The matrix {R} may be {NULL}, in which case only the eigenvalues 
      will be computed.  It may also be the same as {A}, in which case
      {A} will be overwritten by the eigenvectors; otherwise {A} is not 
      altered. */

#endif
