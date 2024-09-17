#ifndef hr2_pmap_affine_encode_H
#define hr2_pmap_affine_encode_H

/* Last edited on 2024-09-17 16:29:05 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <hr2.h>
#include <hr2_pmap.h>
#include <sign.h>

/* A projective map {M} is /affine/ if it takes the hither side of
  {\RT^2} to itself while preserving parallelism between pairs of line
  segments. Its direct matrix {A=M.dir} satisfies {A00 > 0},
  {A10==A20==0}, and {det = A11*A22 - A12*A21 != 0}.
  
  If the determinant {A11*A22-A12*A21} is positive, {M} preserves orientation
  of every point triple, otherwise it reverses it.
  
  In either case, the affine map takes the hither origin to a point with
  Cartesian coordinates {(A01,A02)/A00} (the /displacement/ or {M}), and
  the point difference vectors {(1,0)} {(0,1)} to the point difference
  vectors {(A11,A12)/A00} and {(A21,A22)/A00}, respectively, which are
  not collinear. */

void hr2_pmap_affine_encode(hr2_pmap_t *M, double y[]);
  /*  Converts an affine or neg-affine map {M} to six parameters
    {y[0..5]}. Parameters {y[0..1]} are the Cartesian coordinates of the
    image of the hither point {(0,0)}. Parameters {y[2..3]} are the
    Cartesian coordinates of the image {u} of the difference vector
    {(1,1)}. Parameters {y[4..5]} are a non-linear encoding the image
    {v} of the difference vector {(-1,1)} that is insensitive to
    handedness of the map.
    
    Specifically, let {va,vb} be the components of {v} parallel and
    orthogonal to {u}. Parameter {y[4]} is {va}. Parameter {y[5]} is
    {|vb| - 1/|vb|}. Thus the handedness of {M} (which is the sign of
    {vb}) is lost. */

void hr2_pmap_affine_decode(double y[], hr2_pmap_t *M);
  /* Stores into {M} the affine map with parameters {y[0..5]} and
    positive handedness. Namely, the map will take the hither Cartesian
    origin {(0,0)} to hither Cartesian point {(y[0],y[1])}, the
    difference vector {(1,1)} to {(y[2],y[3])}, and
    the difference vector {(-1,1)} to a nonlinear map
    of {y[3..4]} that always gives a map with positive handedness.

    The functions are inverses in the sense that {decode(y',M);
    encode(M,y')} will set {y''} to {y'}, if {D(y')} is positive, and {encode(M,y);
    decode(y,M)} will have no effect on {M}, provided that {M} is an affine 
    map with positive handedness (except possibly for roundoff errors, if
    {A00} is not 1). */

#endif
