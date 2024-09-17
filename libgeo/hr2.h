/* Oriented projective geometry in two dimensions. */
/* Last edited on 2024-09-17 16:18:23 by stolfi */ 
   
#ifndef hr2_H
#define hr2_H

/* Based on HR2.i3, created 1994-05-04 by J. Stolfi. */

#define _GNU_SOURCE

#include <stdint.h>

#include <r2.h>
#include <r2x2.h>
#include <r3.h>
#include <r3x3.h>

#include <sign.h>

typedef struct hr2_point_t { r3_t c; } hr2_point_t; /* {c.c[0..2]} are the points's coordinates {[w,x,y]}. */
typedef struct hr2_line_t { r3_t f; } hr2_line_t;   /* {f.c[0..2]} are the line's coefficients {<W,X,Y>}. */
  
hr2_point_t hr2_from_r2(r2_t *c);
  /* Point on the ``hither'' half of the two-sided plane (i.e, with positive weight)
    whose Cartesian coordinates are {c}. */
    
r2_t r2_from_hr2(hr2_point_t *p);
  /* Cartesian coordinates of point {p} (which must be finite). */

hr2_point_t hr2_point_at_infinity(r2_t *dir);
  /* The point at infinity whose direction, as seen from any hither
    point, is {*dir}. */

double hr2_pt_pt_diff(hr2_point_t *p, hr2_point_t *q);
  /* Distance between {p} and {q} in the spherical model; that is,
     angle between the vectors {p.c} and {q.c} in {\RR^3}, in radians. */

sign_t hr2_side(hr2_point_t *p, hr2_line_t *L); 
  /* Returns sign of point {p} relative to line {L}: 
    0 on the line, +1 in positive halfplane, -1 in negative halfplane.
    May give inconsistent results for points very close to the line. */

sign_t hr2_orient(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r);
  /* Returns the orientation (turning sense) of the triangle {p q r}.  
    Specifically, the triangle {[1 0 0] [0 1 0] [0 0 1]} 
    has positive orientation. Returns 0 iff the points are collinear.
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
     
double hr2_dist(hr2_point_t *a, hr2_point_t *b);
  /* Euclidean distance between {a} and {b}, which must have weights of the same sign. */
    
double hr2_dist_sqr(hr2_point_t *a, hr2_point_t *b);
  /* Euclidean distance squared between {a} and {b}, which must weights of the same sign. */
     
hr2_point_t hr2_point_mix(double pt, hr2_point_t *p, double qt, hr2_point_t *q);
  /* The point whose homogeneous coordinates are the linear combination
    of those of {p} and {q}, with coefficients {pt} and {qt},
    respectively. */

void hr2_L_inf_normalize_point(hr2_point_t *p);
void hr2_L_inf_normalize_line(hr2_line_t *L);
  /* Scales the homogeneous coordinates/coefficients
    of the given object by a positive factor so that the
    maximum absolute value among them is 1.  If all coordinates
    are zero, they are turned into NaNs. */

hr2_point_t hr2_point_throw(void);
  /* Returns a point with homogeneous coordinates uniformly
    distributed over the unit sphere of {\RR^3}. */
   
hr2_line_t hr2_line_throw(void);
  /* Returns a line with homogeneous coefficienys uniformly
    distributed over the unit sphere of {\RR^3}. */
 

#endif
