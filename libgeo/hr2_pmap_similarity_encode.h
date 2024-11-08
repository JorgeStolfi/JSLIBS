#ifndef hr2_pmap_similarity_encode_H
#define hr2_pmap_similarity_encode_H

/* Last edited on 2024-11-03 06:14:59 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr2.h>
#include <hr2_pmap.h>
#include <sign.h>

/* A /similarity/ is a projective map {M} that takes the hither side of
  {\RT^2} to itself while preserving angles between all triples of finite points.
  Its direct matrix {A=M.dir} satisfies {A00 > 0}, {A10==A20==0},
  {A11^2+A12^2 == A21^2+A22^2}, and {A11*A21+A21*A22 == 0}.
  
  If the determinant factor {A11*A22-A12*A21} is positive, {M} can be
  decomposed into a combination of a uniform scaling {S} by a positive
  factor, a rotation {R}, and a translation {T}; otherwise it is the
  same thing preceded by a swap of the Cartesian {X} and {Y} coordinates.
  
  In either case, the similarity takes the hither origin to a point with
  Cartesian coordinates {(A01,A02)/A00}, and the point difference
  vectors {(1,0)} {(0,1)} to the point difference vectors
  {(A11,A12)/A00} and {(A21,A22)/A00}, respectively; which are
  orthogonal and have equal nonzero length.
  
  The space {\RN_{+}} of all positive similarities is homeomorphic to
  {\RR\times\RS^1\times\RR^2}, where the first {\RR} factor accounts for
  the scaling, the second factor {\RS^1} for the rotation, and the third
  factor {\RR^2} for the translation. */

void hr2_pmap_similarity_encode(hr2_pmap_t *M, double y[]);
  /* Assumes that {M} is positive similarity map, and converts it to four
    parameters {enc(M) = y[0..3]}, which are the scaling factor the angle of rotation,
    and the displacement vector of the translation.
  
    Specifically, the procedure sets sets {y[0]} to {0.5*(f-1/f)}, where
    {f} is the (positive) scaling factor, namely {hypot(A11,A12)/A00}
    which should be the same as {hypot(A21,A22)/A00}. Then it sets
    {y[1]} to the CCW angle (in radians) between the {X}-axis direction
    and the image of the difference vector {(1,0)}, that is,
    {(A11,A12)/A00}. Finally, sets {y[2]} and {y[3]} to the Cartesian
    coordinates of the image of the origin {(0,0)}, namely {A01/A00} and
    {A02/A00}.
    
    The rotation angle {y[1]} is always in the range {[-PI _ +PI)}. */

void hr2_pmap_similarity_decode(double y[], hr2_pmap_t *M); 
  /* Stores into {M} the positive similarity map {dec(y)} with
    parameters {y[0..3]}. 
    
    The resulting map will scale the Cartesian {X} and {Y} by 
    {f=y[0]+hypot(y[0],1)}, rotate the vector {(1,0)} by {y[1]} radians CCW, and
    take the hither origin {(0,0)} to the hither point {(y[2],y[3])},.

    The function is the inverse of {encode} in the sense that
    {encode(M,y); decode(y,M)} will have no net effect on {M} if {M} is
    a positive similarity; and {decode(y,M); encode(M,y)} will have no
    net effect on {y}, except that {y[1]} will be reduced to the range
    {(-PI _ +PI]} module {2*PI}. However, allowance must be made in both
    for possible floating point roundoff errors. */

#endif
