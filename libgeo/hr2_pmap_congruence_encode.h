#ifndef hr2_pmap_congruence_encode_H
#define hr2_pmap_congruence_encode_H

/* Last edited on 2024-12-05 10:26:59 by stolfi */

#include <stdint.h>

#include <hr2.h>
#include <hr2_pmap.h>
#include <sign.h>

/* A /congruence/ is a projective map {M} that takes the hither side of
  {\RT^2} to itself while preserving distances between all pairs of
  finite points. Its direct matrix {A=M.dir} satisfies {A00 > 0},
  {A10==A20==0}, {A11^2+A12^2 == A21^2+A22^2 == A00^2}, and {A11*A21 +
  A21*A22 == 0}.
  
  If the determinant {A11*A22-A12*A21} is positive, {M} is a combination
  of a rotation {R} and a translation {T}; otherwise it is a combination
  of a reflection {L} about the {X}-axis, a rotation {R}, and a
  translation {T}.
  
  In either case, the congruence takes the hither origin to a point with
  Cartesian coordinates {(A01,A02)/A00}, and the point difference
  vectors {(1,0)} {(0,1)} to the point difference vectors
  {(A11,A12)/A00} and {(A21,A22)/A00}, respectively; which are
  orthogonal and have unit length.
  
  The space {\RN_{+}} of all positive congruences is homeomorphic to
  {\RS^1\times \RR^2}, where the first factor {\RS^1} accounts for the
  rotation and the second factor {\RR^2} for the translation. */

void hr2_pmap_congruence_encode(hr2_pmap_t *M, double y[]);
  /* Assumes that {M} is a positive congruence map, and converts it to
    three parameters {enc(M) = y[0..2]}, which are the angle of
    rotation and the displacement vector.
  
    Specifically, the procedure sets {y[0]} to the CCW angle (in radians)
    from the {X} axis to the image of the difference vector {(1,0)},
    that is, {(A11,A12)/A00}. Then it sets {y[1]} and {y[2]} to the Cartesian
    coordinates {A01/A00} and {A02/A00} of the image of the hither
    origin {(0,0)}.
    
    The rotation angle {y[0]} is always in the range {[-PI _ +PI)}. */

void hr2_pmap_congruence_decode(double y[], hr2_pmap_t *M);
  /* Stores into {M} the positive congruence map {dec(y)} with
    parameters {y[0..2]}. The resulting map will rotate the difference
    vectors {(1,0)} and {(0,1)} by {y[0]} radians CCW, and take the hither origin
    {(0,0)} to the hither point {(y[1],y[2])}.
    
    The function is the inverse of {encode} in the sense that
    {encode(M,y);decode(y,M)} will have no net effect on {M} if {M} is a
    positive congruence, and {decode(y,M);encode(M,y)} will have no net
    effect on {y}, except that {y[0]} will be reduced to the range {(-PI
    _ +PI]} module {2*PI}. However, allowance must be made in both for
    possible floating point roundoff errors. */

#endif
