/* Oriented projective geometry in two dimensions. */
/* Last edited on 2020-10-12 14:12:14 by jstolfi */ 
   
#ifndef hr2_H
#define hr2_H

#include <r2.h>
#include <r3.h>
#include <r2x2.h>
#include <r3x3.h>
#include <sign.h>

/* Based on HR2.i3, created 1994-05-04 by J. Stolfi. */

typedef struct hr2_point_t { r3_t c; } hr2_point_t; /* {c.c[0..2]} are the points's coordinates. */
typedef struct hr2_line_t { r3_t f; } hr2_line_t;   /* {f.c[0..2]} are the line's coefficients. */
  
hr2_point_t hr2_from_r2(r2_t *c);
  /* Point on ``hither'' half of the plane (i.e, with positive weight)
    whose Cartesian coordinates are {c}. */
    
r2_t r2_from_hr2(hr2_point_t *p);
  /* Cartesian coordinates of point {p} (which must be finite). */

double hr2_pt_pt_diff(hr2_point_t *p, hr2_point_t *q);
  /* Distance between {p} and {q} in the spherical model; that is,
     angle between the vectors {p.c} and {q.c} in {R^3}, in radians. */

sign_t hr2_side(hr2_point_t *p, hr2_line_t *L); 
  /* Returns sign of point {p} relative to line {L}: 0 on the line,
    +1 in positive halfplane, -1 in negative halfplane.
    May give inconsistent results for points very close to the line. */

sign_t hr2_orient(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r);
  /* Returns the orientation (turning sense) of the triangle {p q r}.  
    Specifically, the triangle {[1 0 0] [0 1 0] [0 0 1]} has positive
    orientation. Returns 0 iff the points are collinear.
    May give inconsistent results for points very close to collinear. */

hr2_line_t hr2_join(hr2_point_t *p, hr2_point_t *q);
  /* Return the line through {p} and {q}. */

hr2_point_t hr2_meet(hr2_line_t *K, hr2_line_t *L);
  /* Return the point common to {K} and {L}. */

r2_t hr2_point_point_dir(hr2_point_t *frm, hr2_point_t *tto);
  /* Direction (a unit-length vector) of point {tto} seen from point {frm}.
    Works even if one of them is at infinity.  Does not work if both
    are at infinity, or coincident, or antipodal. */
    
r2_t hr2_line_dir(hr2_line_t *L);
  /* The direction of line {L}: a unit vector which, on the hither side
    of the plane, is parallel to {L} and travels around its positive
    side in the counterclockwise sense. Assumes {L} is not at
    infinity. */
    
r2_t hr2_line_normal(hr2_line_t *L);
  /* The normal direction of line {L}: a unit vector which,
    on the hither side of the plane, points from {L} into 
    {L}'s positive halfplane.  Assumes {L} is not at infinity. */
    
/* PROJECTIVE MAPS */

typedef struct hr2_pmap_t { r3x3_t dir; r3x3_t inv; } hr2_pmap_t;
  /* A projective map. Field {dir} is the map's matrix, {inv} is its inverse. */

bool_t hr2_pmap_is_identity(hr2_pmap_t *m);
  /* TRUE iff {m} is the identity map (apart from homogeneous scaling). */

hr2_point_t hr2_pmap_point(hr2_point_t *p, hr2_pmap_t *m);
  /* Applies projective map {m} to point {p}. */

hr2_line_t hr2_pmap_line(hr2_line_t *L, hr2_pmap_t *m);
  /* Applies projective map {m} to line {L}. */

hr2_pmap_t hr2_pmap_translation(r2_t *vec);
  /* Returns a projective map that performs a translation by the Cartesian vector {vec}.
    The map is actually affine (has {[1 0 0]} as the first column). */

hr2_pmap_t hr2_pmap_rotation(double ang);
  /* Returns a projective map that performs a rotaton by {ang} radians about the origin.
    The map is actually linear (has {[1 0 0]} as the first column and the first row). */

hr2_pmap_t hr2_pmap_from_mat_and_disp(r2x2_t *mat, r2_t *disp);
  /* Returns an affine map that performs the linear map 
    described by the matrix {mat} followed by translation by {disp}.
    That is, maps the point with Cartesian coordinates {p \in \RR^2} 
    to {p*mat + disp}. */

hr2_pmap_t hr2_pmap_comp(hr2_pmap_t *m, hr2_pmap_t *n);
  /* Returns the composition of {m} and {n}, applied in that order. */

hr2_pmap_t hr2_pmap_inv(hr2_pmap_t *m);
  /* Returns the inverse of map {m}. */

hr2_pmap_t hr2_pmap_from_points(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r, hr2_point_t *u);
  /* Returns a projective map that takes the cardinal points {[1,0,0]},
    {[0,1,0]}, and {[0,0,1]} to {p}, {q}, and {r}, respectively; and
    also some point {s} of the form {[±1,±1,±1]} to {u}. 
    
    The point {s} is unique, and is the signature of {u} relative to
    the ordered triple {p,q,r}.
    
    The procedure fails if the set {p,q,r,u} contains three collinear
    points.  */
    
#endif
