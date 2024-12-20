/* sym_eigen_tred1.h -- tridiagonalization of a symmetric matrix. */
/* Last edited on 2024-12-05 11:16:32 by stolfi */

#ifndef sym_eigen_tred1_H
#define sym_eigen_tred1_H

#include <stdint.h>

void sym_eigen_tred1
  ( uint32_t n,
    double A[],
    double d[],
    double e[],
    double e2[]
  );
  /* This procedure is based on the EISPACK routines "tred1.f" and
    "tred2.f" by Burton S. Garbow, Argonne National Laboratory
    (aug/1983); originally from the ALGOL procedures {tred1} and
    {tred2}, Num. Math. 11, 181-195(1968) by Martin, Reinsch, and
    Wilkinson. See Handbook for Auto. Comp., vol.II-Linear Algebra,
    212-226(1971). Re-implemented in C by Jorge Stolfi, Unicamp
    (dec/2002 and dec/2024).

    This proceduer reduces a real symmetric matrix to a symmetric
    tridiagonal matrix using orthogonal similarity transformations.

    On input

       {n} is the order of the matrix.

       {A} contains the real symmetric input matrix.  only the
         lower triangle of the matrix need be supplied.
         The matrix should be stored by rows; that is, {A[i,j]}
         is assumed to actually be {A[i*n + j]}.

    On output

       {A} contains information about the orthogonal trans-
         formations used in the reduction in its strict lower
         triangle.  the full upper triangle of {A} is unaltered.

       {d} contains the diagonal elements of the tridiagonal matrix.

       {e} contains the subdiagonal elements of the tridiagonal
         matrix in its last n-1 positions.  {e[0]} is set to zero.

       {e2} contains the squares of the corresponding elements of {e}.
         {e2} may coincide with {e} if the squares are not needed. */

#endif
