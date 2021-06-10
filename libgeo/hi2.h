/* Oriented projective geometry in two dimensions with integer homogeneous coordinates. */
/* Last edited on 2021-06-09 20:29:07 by jstolfi */ 
   
#ifndef hi2_H
#define hi2_H

#define _GNU_SOURCE
#include <stdint.h>

#include <i2.h>
#include <i3.h>
/* #include <i3x3.h> */
#include <sign.h>
#include <urat64.h>

typedef struct hi2_point_t { i3_t c; } hi2_point_t; /* {c.c[0..2]} are the points's coordinates. */
typedef struct hi2_line_t { i3_t f; } hi2_line_t;   /* {f.c[0..2]} are the line's coefficients. */

hi2_point_t hi2_from_i2(i2_t *c);
  /* Point on ``hither'' half of the plane (i.e, with positive weight)
    whose Cartesian coordinates are {c}.
    Input MMax: {2^31-1}. Output range: {-M .. +M}. */
    
sign_t hi2_side(hi2_point_t *p, hi2_line_t *L); 
  /* Returns sign of point {p} relative to line {L}: 0 on the line,
    +1 in positive halfplane, -1 in negative halfplane.
    Input MMax: {1753413056 > 1672*2^20}. */

sign_t hi2_orient(hi2_point_t *p, hi2_point_t *q, hi2_point_t *r);
  /* Returns the orientation (turning sense) of the triangle {p q r}.  
    Specifically, the triangle {[1 0 0] [0 1 0] [0 0 1]} has positive
    orientation. Returns 0 iff the points are collinear.
    Input MMax: {1518500249 > 1448*2^20}. */

hi2_line_t hi2_join(hi2_point_t *p, hi2_point_t *q);
  /* Return the line through {p} and {q}. 
    Input MMax: {2^15-1 = 32767}. Output range: {-2*M^2 .. +2*M^2}. */

hi2_point_t hi2_meet(hi2_line_t *K, hi2_line_t *L);
  /* Return the point common to {K} and {L}. 
    Input MMax: {2^15-1 = 32767}. Output range: {-2*M^2 .. +2*M^2}. */

i2_t hi2_point_point_dir(hi2_point_t *frm, hi2_point_t *tto);
  /* Returns a vector pointing along the shortest path from point {frm}
    to point {tto}.  Works even if one of them is at infinity.  
    Does not work if both are at infinity, or coincident, or antipodal.
    Input MMax: {2^15-1 = 32767}. Output range: {-2*M^2 .. +2*M^2}. */
    
i2_t hi2_line_dir(hi2_line_t *L);
  /* Returns a vector which, on the hither side
    of the plane, is parallel to {L} and travels around its positive
    side in the counterclockwise sense. Assumes {L} is not at infinity.
    Input MMax: {2^31-1}. Output range: {-M .. +M}. */
    
i2_t hi2_line_normal(hi2_line_t *L);
  /* The normal direction of line {L}: a unit vector which,
    on the hither side of the plane, points from {L} into 
    {L}'s positive halfplane.  Assumes {L} is not at infinity.
    Input MMax: {2^31-1}. Output range: {-M .. +M}. */
    
urat64_t hi2_dist_sqr(hi2_point_t *a, hi2_point_t *b);
  /* The Euclidean distance squared from {a} to {b} or to the 
    antipode of {b}, whichever is shorter. 
    Input MMax: {38967 > 38*2^10}. 
    Output range: num {0 .. 8*M^4}, den {0 .. M^4}. */

sign_t hi2_in_circle(hi2_point_t *a, hi2_point_t *b, hi2_point_t *c, hi2_point_t *d);
  /* Returns {-1,0,+1} depending on whether {d} is inside the circle  
    through the points vectors {a,b,c}. 
    Input MMax: {2^15-1 = 32767}. */


// /* PROJECTIVE MAPS */
// 
// typedef struct hi2_pmap_t { i3x3_t dir; i3x3_t inv; } hi2_pmap_t;
//   /* A projective map. Field {dir} is the map's matrix, {inv} is its inverse. */
// 
// bool_t hi2_pmap_is_identity(hi2_pmap_t *m);
//   /* TRUE iff {m} is the identity map (apart from homogeneous scaling). */
// 
// hi2_point_t hi2_map_point(hi2_point_t *p, hi2_pmap_t *m);
//   /* Applies projective map {m} to point {p}. */
// 
// hi2_line_t hi2_map_line(hi2_line_t *L, hi2_pmap_t *m);
//   /* Applies projective map {m} to line {L}. */
// 
// hi2_pmap_t hi2_comp_map(hi2_pmap_t *m, hi2_pmap_t *n);
//   /* Returns the composition of {m} and {n}, applied in that order. */
// 
// hi2_pmap_t hi2_inv_map(hi2_pmap_t *m);
//   /* Returns the inverse of map {m}. */
// 
// hi2_pmap_t hi2_pmap_from_points(hi2_point_t *p, hi2_point_t *q, hi2_point_t *r, hi2_point_t *u);
//   /* Returns a projective map that takes the cardinal points {[1,0,0]},
//     {[0,1,0]}, and {[0,0,1]} to {p}, {q}, and {r}, respectively; and
//     also some point {s} of the form {[±1,±1,±1]} to {u}. 
//     
//     The point {s} is unique, and is the signature of {u} relative to
//     the ordered triple {p,q,r}.
//     
//     The procedure fails if the set {p,q,r,u} contains three collinear
//     points.  */
    
#endif
