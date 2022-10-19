/* rmxn_transform_quadratic.h --- restricts a quadratic form to a linear subspace */
/* Last edited on 2021-12-24 20:31:50 by stolfi */

#ifndef rmxn_transform_quadratic_H
#define rmxn_transform_quadratic_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <rn.h>
#include <rmxn.h>

void rmxn_transform_quadratic(int32_t n, double E[], double e[], int32_t m, double U[], double F[], double f[]);
  /* The parameters {E}, {U}, and {F} must be matrices of size {n×n}, {m×n}, and {m×m},
    respectively, linearized by rows. 
    
    The procedure interprets {E} and {e} as defining a quadratic form
    {\QE} on {\RR^n} that takes any row {n}-vector {x} to {\QE(x) =
    x*E*diag(e)*E'*x'}. It also interprets {U} as the matrix of a linear
    map {\LU} from {\RR^m} to {\RR^n} such that {\LU(z) = z*U} for any row
    {m}-vector {z}.
    
    The procedure stores in {F} and {f} the orthonormal {m×m} matrix and
    the {m}-vector such that {\QF(z) = z*F*diag(f)*F'*z'} for any row
    {m}-vector {z}, where {\QF} is the quadratic form defined on {\RR^m}
    such that {\QF(z) = \QE(\LU(z))}. That is, {F} and {f} will be the
    eigenvectors and eigenvalues of {U*E*diag(e)*E'*U'}. The vector {f}
    will be sorted in ascending order.
    
    Said another way, let {\EE} ne the ellipsoid of {\RR^n} with radius
    {\rE[i]} along axis {E[i]}, and {e[i]} be {1/\rE[i]^2} for all {i}.
    Then {sqrt(\QE)} is the norm of {\RR^n} whose unit ball is {\EE}.
    Further, let {\EU} be the intersection of {\EE} with the linear
    subspace {\SU = \LU(\RR^m)}, that is, with the row span of {U}. In
    that case, {sqrt(\QF)} will be the norm of {\RR^m} whose unit ball
    is the ellipsoid {\EF}, such that {\EU = \LU(\EF)}; and {f[k]}, for
    {k} in {0..m-1}, will be {1/\rF[k]^2} where {\rF[k]} is the radius
    of {\EF} along axis {F[k]}. */

#endif

