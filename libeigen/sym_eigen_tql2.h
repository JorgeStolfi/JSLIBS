/* sym_eigen_tql2.h -- eigenvalues and eigenvectors of a tridiagonal matrix */
/* Last edited on 2024-12-05 11:11:37 by stolfi */

#ifndef sym_eigen_tql2_H
#define sym_eigen_tql2_H

#include <stdint.h>

void sym_eigen_tql2
  ( uint32_t n,
    double d[],
    double e[],
    double z[],
    uint32_t *nev_P
  );
  /* Finds the right eigenvalues and eigenvectors of a symmetric
    tridiagonal matrix {T} by the QL method. 
    
    If the matrix {T} vas derived from a full symmetric matrix {A} by
    an orthogonal similarity map {z} (see {sym_eigen_tridiagonalize}), this
    procedure will find the eigenvectors of {A}.

    On input,

       {n} is the size (rows and columns) of the matrices {T} and {z}.

       {d} must contain the diagonal elements of the input matrix {T}.

       {e} must contain the subdiagonal elements of the input matrix {T}
         in its last {n-1} positions (i.e. {e[i] == T[i,i-1]},
         for {i=1..n-1}).  The value of {e[0]} is ignored.

       {z} must contain the orthogonal similarity matrix which
         transformed the original matrix {A} into {T}, i.e. the matrix
         such that {z*A*(z^t) == T}. (To obtain the eigenvectors of
         {T} itself, {z} must be set to the identity matrix.)

     On output,

       {*nev_P} is the number {nev} of eigenvalues and eigenvectors that
         could be computed, which is {n} in successful run. If less than
         {n}, it means that the computation of {d[nev]} failed to
         converge in 30 iterations, thus {d[nev..n-1]} are undefined.

       {d[0..*nev-1]} contains the eigenvalues that could be 
         determined, in ascending order of signed value.
         If {nev < n}, they will be some subset of the eigenvalues,
         but not necessarily the {nev} smallest ones in signed value.

       {e} has been destroyed.

       {z[0..*nev-1,0..n-1]} are orthonormal eigenvectors of the matrix
         {A}, corresponding to the eigenvalues {d[0..*nev-1]}, stored as
         rows. If {S} is {z} truncated to its first {nev} rows, then
         {S*A*(S^t) == diag(d[0..nev-1])}. */

#endif
