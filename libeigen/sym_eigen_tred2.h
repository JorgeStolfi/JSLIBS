/* sym_eigen_tred2.h -- tridiagonalization of a symmetric matrix. */
/* Last edited on 2024-12-05 11:16:40 by stolfi */

#ifndef sym_eigen_tred2_H
#define sym_eigen_tred2_H

#include <stdint.h>

void sym_eigen_tred2
  ( uint32_t n,
    double A[],
    double d[],
    double e[],
    double Z[]
  );
  /* This procedure is based on the EISPACK routines "tred1.f" and
    "tred2.f" by Burton S. Garbow, Argonne National Laboratory
    (aug/1983); originally from the ALGOL procedures {tred1} and
    {tred2}, Num. Math. 11, 181-195(1968) by Martin, Reinsch, and
    Wilkinson. See Handbook for Auto. Comp., vol.II-Linear Algebra,
    212-226(1971). Re-implemented in C by Jorge Stolfi, Unicamp
    (dec/2002 and dec/2024). 

    This procedure reduces a real symmetric matrix to a symmetric
    tridiagonal matrix using and accumulating orthogonal similarity
    transformations.

    On input

       {n} is the order of the matrix.

       {A} contains the real symmetric input matrix.  only the
         lower triangle of the matrix need be supplied.
         The matrix should be stored by rows; that is, {A[i,j]}
         is assumed to actually be {A[i*n + j]}.

    On output

       {d} contains the diagonal elements of the tridiagonal matrix.

       {e} contains the subdiagonal elements of the tridiagonal
         matrix in its last n-1 positions.  {e[0]} is set to zero.

       {Z} contains the orthogonal transformation matrix
         produced in the reduction.

       {A} and {Z} may coincide.  if distinct, {A} is unaltered. */

#endif
