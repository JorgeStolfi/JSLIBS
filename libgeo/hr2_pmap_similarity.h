/* Similarity maps in 2D oriented projective plane {\RT^2}. */
/* Last edited on 2024-12-05 10:27:03 by stolfi */ 
   
#ifndef hr2_pmap_similarity_H
#define hr2_pmap_similarity_H

#include <stdint.h>

#include <sign.h>
#include <bool.h>
#include <r2.h>
#include <hr2.h>
#include <hr2_pmap.h>

hr2_pmap_t hr2_pmap_similarity_from_two_points(r2_t *p, r2_t *q, sign_t sgn);
  /* Returns a projective map that is a Cartesian similarity taking the
    points {[1,0,0]}, and {[1,1,0]} to {p} and {q}, respectively. The
    points must be distinct.
    
    A similarity (or Euclidean transformation) is a map from {\RR^2} to
    {\RR^2} that preserves ratios between distances, and therefore
    preserves all angles. It is a special case of affine map. It
    preserves the sign of homogeneous coordinate 0 (weight).   Note that a
    uniform scaling by a negative factor is equivalent to a rotation by
    180 degrees.
    
    The {sgn} parameter should be nonzero. If {sgn} is {+1}, the map
    will preserve handedness, i. e. will be a combination of a rotation
    plus a translation. If {sgn} is {-1}, the map will reverse
    handedness; it will be a combination of rotation, translation, and
    mirroring about a line.  */

#endif
