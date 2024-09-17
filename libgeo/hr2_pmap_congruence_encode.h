#ifndef hr2_pmap_congruence_encode_H
#define hr2_pmap_congruence_encode_H

/* Last edited on 2024-09-17 16:29:23 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr2.h>
#include <hr2_pmap.h>
#include <sign.h>

/* A /congruence/ is a projective map {M} that takes the hither side of
  {\RT^2} to itself while preserving distances between all pairs of finite points.
  Its direct matrix {A=M.dir} satisfies {A00 > 0}, {A10==A20==0},
  {A11^2+A12^2 == A21^2+A22^2 == A00^2}, and {A11*A21 + A21*A22 == 0}.
  
  If the determinant {A11*A22-A12*A21} is positive, {M} is a combination
  of a rotation {R} and a translation {T}; otherwise it is a combination of 
  a reflection {L} about the {X}-axis, a rotation {R}, and a translation {T}.
  
  In either case, the congruence takes the hither origin to 
  a point with Cartesian coordinates {(A01,A02)/A00}, and 
  the point difference vectors {(1,0)} {(0,1)} to the point difference vectors
  {(A11,A12)/A00} and {(A21,A22)/A00}, respectively;
  which are orthogonal and have unit length. */

void hr2_pmap_congruence_encode(hr2_pmap_t *M, double y[]);
  /* Assumes that {M} is a congruence map, and converts it to three
    parameters {y[0..2]}, which are the Cartesian coordinates of the image of the
    origin {(0,0)} and the angle of rotation of the vector {(1,0)}.
  
    Specifically, the procedure sets {y[0]} and {y[1]} to the Cartesian
    coordinates {A01/A00} and {A02/A00} of the image of the hither
    origin {[1,0,0]}. Then it sets {y[2]} to the CCW angle (in radians)
    from the {X} axis to the image of the vector {(1,0)} or {(0,1)}, namely
    to {(A11,A12)} or {(A21,A22)}, depending on the handedness of {M}.
    
    Note that the handedness of {M} is lost. */

void hr2_pmap_congruence_decode(double y[], hr2_pmap_t *M);
  /* Stores into {M} a congruence map with parameters {y[0..2]}
    and positive handedness.  The resulting map will take the hither origin {(0,0)}
    to the hither point {(y[0],y[1])}, and rotate the difference vectors {(1,0)} and {(0,1)}
    by {y[2]} radians CCW.
    
    The functions are inverses in the sense that {decode(y',M);
    encode(M,y')} will set {y''} to {y'}, and {encode(M,y);
    decode(y,M)} will have no effect on {M}, if {M} has positive
    handedness (except possibly for roundoff errors, if {A00} is not 1). */

#endif
