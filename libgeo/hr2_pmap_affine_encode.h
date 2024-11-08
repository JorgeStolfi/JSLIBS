#ifndef hr2_pmap_affine_encode_H
#define hr2_pmap_affine_encode_H

/* Last edited on 2024-11-02 23:40:39 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr2.h>
#include <hr2_pmap.h>
#include <sign.h>

/* A projective map {M} is /affine/ if it takes the hither side of
  {\RT^2} to itself while preserving parallelism between pairs of line
  segments. Its direct matrix {A=M.dir} satisfies {A00 > 0},
  {A10==A20==0}, and {A11*A22-A12*A21 != 0}. 
  If the determinant factor {A11*A22-A12*A21} is positive, {M} preserves
  orientation of every point triple, otherwise it reverses it.
  
  A positive affine map {M} can be decomposed into a combination of a
  shearing and {Y}-scaling map {H} that keeps the points {(0,0)} and
  {(1,0)} fixed while mapping {(0,1)} to some vector {(vX,vY)} with
  {vY>0}; a uniform scaling {S} by a positive factor; a rotation {R},
  and a translation {T}. If the map is negative, it is the same thing,
  preceded by a swap of the {X} and {Y} coordinates.
  
  In any case, an affine map takes the hither origin to a point with
  Cartesian coordinates {(A01,A02)/A00} (the /displacement/ or {M}), and
  the point difference vectors {(1,0)} {(0,1)} to the point difference
  vectors {(A11,A12)/A00} and {(A21,A22)/A00}, respectively, which are
  not collinear.
  
  The space {\RN_{+}} of all positive affine maps is homeomorphic to
  {\RR^2\times\RR\times\RS^1\times\RR^2}, where the first factor {\RR^2}
  accounts for the shearing, the second one {\RR} for the scaling, the
  third one {\RS^1} for the rotation, and the last one {\RR^2} for the
  translation. */

void hr2_pmap_affine_encode(hr2_pmap_t *M, double y[]);
  /* Assumes that {M} is a positive affine map, and converts it to six
    parameters {enc(M) = y[0..5]}, which define the shearing, scaling,
    rotation, and translation component of the map.
    
    Specifically, let {u=(A11,A12)/A00} and {v=(A21,A22)/A00} be the
    images of the difference vectors {(1,0)} and {(0,1)}, respectively.
    Let {u'} be the vector {u} rotated 90 degrees CCW. Let {c} and {s}
    be coefficients such that {v = c*u + s*u'}. Note that {c} is
    arbitrary but {s} must be positive. The procedure sets {y[0]} and
    {y[1]} to {c} and {0.5*(s-1/s)}. Then it sets {y[2]} to
    {0.5*(f-1/f)} where {f} is the length of {u}, and {y[3]} to the CCW
    angle in radians between {(1,0)} and {u}. Finally it sets {y[4]} and
    {y[5]} to the the Cartesian coordinates of the image of the hither
    point {(0,0)}.
    
    The rotation angle {y[3]} is always in the range {[-PI _ +PI)}. */

void hr2_pmap_affine_decode(double y[], hr2_pmap_t *M);
  /* Stores into {M} the positive affine map {dec(y)} with parameters {y[0..5]}.
    The resulting map will first apply a shearing and {Y}-scaling 
    map that keeps points {(0,0)} and {(1,0)} fixed while
    taking the direction vector {(0,1)} to the point 
    {(c,s)} where {c} is {y[0]} and {s} is {y[1]+hypot(y[1],1)}.
    Then it applies a uniform scaling by a factor
    {f=y[2]+hypot(y[2],1)} and a rotation by a CCW angle of {y[3]} radians
    about the origin.  Finally it applies a translation that takes the hither Cartesian
    origin {(0,0)} to hither Cartesian point {(y[0],y[1])}.

    The functionis the inverse of {encode} in the sense that
    {encode(M,y); decode(y,M)} will have no net effect on {M} if {M} is
    a positive affine map; and {decode(y,M); encode(M,y)} will have no
    effect on {y}, except that {y[3]} will be reduced to the range {(-PI
    _ +PI]} module {2*PI}. However, allowance must be made in both for
    possible floating point roundoff errors. */

#endif
