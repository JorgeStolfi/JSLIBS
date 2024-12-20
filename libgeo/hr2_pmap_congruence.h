/* Congruence maps in 2D oriented projective plane {\RT^2}. */
/* Last edited on 2024-12-05 10:26:56 by stolfi */ 
   
#ifndef hr2_pmap_congruence_H
#define hr2_pmap_congruence_H

#include <stdint.h>

#include <sign.h>
#include <r2.h>
#include <hr2_pmap.h>

hr2_pmap_t hr2_pmap_congruence_from_point_and_dir(r2_t *p, r2_t *u, sign_t sgn);
  /* Returns a projective map that is actually a Cartesian congruence
    taking the origin {[1,0,0]} to the Cartesian point {p}
    the direction of the {X}-axis to the direction vector {u}.
    The length of {u} is ignored, but must not be zero.
    
    A Cartesian congruence (or isometry) is a map of {\RR^2} to {\RR^2}
    that preserves all distances, and therefore also all angles. It is a
    special case of similarity and affine map. It preserves the sign
    of homogeneous coordinate 0 (weight). 
    
    The {sgn} parameter should be nonzero. If {sgn} is {+1}, the map
    will preserve handedness, i. e. will be a combination of a rotation
    plus a translation. If {sgn} is {-1}, the map will reverse
    handedness; it will be a combination of rotation, translation, and
    mirroring about a line. */

#endif
