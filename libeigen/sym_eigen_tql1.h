/* sym_eigen_tql1.h -- eigenvalues and eigenvectors of a tridiagonal matrix */
/* Last edited on 2024-12-05 11:17:16 by stolfi */

#ifndef sym_eigen_tql1_H
#define sym_eigen_tql1_H

#include <stdint.h>

void sym_eigen_tql1
  ( uint32_t n,
    double d[],
    double e[],
    uint32_t *nev_P
  );
  /* Finds the right eigenvalues and eigenvectors of a symmetric
    tridiagonal matrix {T} by the QL method. 

    On input,

       {n} is the size (rows and columns) of the matrix {T}.

       {d[0..n-1]} must be the diagonal elements of the input matrix {T}.

       {e[1..n-1]} must be the the subdiagonal elements of the input matrix {T}
         (i.e. {e[i] == T[i,i-1] == T[i-1,i]} for {i=1..n-1}). 
         The value of {e[0]} is ignored.

     On output,

       {*nev_P} is the number {nev} of eigenvalues that were
         successfully computed, which is {n} after a successful run. if
         less than {n}, it means that the computation of {d[nev]} failed
         to converge in 30 iterations, thus {d[nev..n-1]} are undefined.

       {d[0..*nev-1]} contains the eigenvalues that could be determined,
         in ascending order of signed value. If {nev < n}, they will be
         some subset of the eigenvalues, but not necessarily the {nev}
         smallest ones in signed or absolute value.

       {e} has been destroyed. */

#endif
