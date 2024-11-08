/* Oriented projective geometry in three dimensions. */
/* Last edited on 2024-11-03 07:36:54 by stolfi */ 

#ifndef hr3_H
#define hr3_H

/* Based on HR3.i3, created 1993-04-18 by Marcos C. Carrard. */
/* Also based on H3.pas by J. Stolfi. */
   
#define _GNU_SOURCE

#include <stdint.h>

#include <r3.h>
#include <r3x3.h>
#include <r4.h>
#include <r4x4.h>
#include <r6.h>

#include <sign.h>

typedef struct hr3_point_t { r4_t c; } hr3_point_t; /* {c.c[0..3]} are point's coords {[w,x,y,z]}. */
typedef struct hr3_plane_t { r4_t f; } hr3_plane_t; /* {f.c[0..3]} are plane's coeffs {<W,X,Y,Z>}. */
typedef struct hr3_line_t  { r6_t k; } hr3_line_t;  /* {k.c[0..5]} are Plücker coords {[wx,wy,xy,wz,xz,yz]} */
  
hr3_point_t hr3_from_r3(r3_t *c);
  /* Point on the ``hither'' half of the two-sided space (i.e, with positive weight)
    whose Cartesian coordinates are {c}. */
    
r3_t r3_from_hr3(hr3_point_t *p);
  /* Cartesian coordinates of point {p} (which must be finite). */
  
r3_t r3_from_hr3_nan(hr3_point_t *p);
  /* Converts a point from homogeneous coordinates {p} to Cartesian
    coordinates {c}.  If {p} is at infinity, returns {(NAN,NAN,NAN)}. */

hr3_point_t hr3_point_at_infinity(r3_t *dir);
  /* The point at infinity whose direction, as seen from any hither
    point, is {*dir}. */

double hr3_pt_pt_diff(hr3_point_t *p, hr3_point_t *q);
  /* Distance between {p} and {q} in the spherical model; that is,
     angle between the vectors {p.c} and {q.c} in {\RR^4}, in radians. */

sign_t hr3_side(hr3_point_t *p, hr3_plane_t *Q); 
  /* Returns the position of point {p} relative to plane {Q}: 
    0 on the plane, +1 in positive halfspace, -1 in negative halfspace.
    May give inconsistent results for points very close to the plane. */

sign_t hr3_orient(hr3_point_t *p, hr3_point_t *q, hr3_point_t *r, hr3_point_t *s);
  /* Returns the orientation (handedness) of the tetrahedron {p q r s}.
    Specifically, the tetrahedron {[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 1]}
    has positive orientation. Returns 0 iff the points are coplanar.
    May give inconsistent results for points very close to coplanar. */
    
hr3_line_t hr3_line_from_two_points(hr3_point_t *p, hr3_point_t *q);
  /* The line from {p} to {q}. */

hr3_plane_t hr3_plane_from_three_points(hr3_point_t *p, hr3_point_t *q, hr3_point_t *r);
  /* The plane through {p}, {q}, and {r}. */
    
hr3_plane_t hr3_plane_from_line_and_point(hr3_line_t *L, hr3_point_t *r);
  /* The plane through {L} and {r}. */
    
hr3_line_t hr3_line_from_two_planes(hr3_plane_t *P, hr3_plane_t *Q);
  /* The line where {P} meets {Q}. */

hr3_point_t hr3_point_from_three_planes(hr3_plane_t *P, hr3_plane_t *Q, hr3_plane_t *R);
  /* The point where {P}, {Q}, and {R} meet. */
    
hr3_point_t hr3_point_from_line_and_plane(hr3_line_t *L, hr3_plane_t *R);
  /* The point where {L} meets {R}. */

r3_t hr3_point_point_dir(hr3_point_t *frm, hr3_point_t *tto);
  /* Direction (a unit-length vector) of point {tto} seen from point {frm}.
    Works even if one of them is at infinity.  Does not work if both
    are at infinity, or coincident, or antipodal. */
      
r3_t hr3_plane_normal(hr3_plane_t *P);
  /* The normal direction of plane {P}, on the hither side, pointing into 
    {P}'s positive halfspace.  Assumes {P} is not at infinity. */

double hr3_dist(hr3_point_t *a, hr3_point_t *b);
  /* Euclidean distance between {a} and {b}, which must have weights of the same sign. */
    
double hr3_dist_sqr(hr3_point_t *a, hr3_point_t *b);
  /* Euclidean distance squared between {a} and {b}, which must weights of the same sign. */
    
r3_t hr3_line_dir(hr3_line_t *L);
  /* The direction along line {L}.  Assumes {L} is not at infinity. */
    
hr3_point_t hr3_point_mix(double pt, hr3_point_t *p, double qt, hr3_point_t *q);
  /* The point whose homogeneous coordinates are the linear combination
    of those of {p} and {q}, with coefficients {pt} and {qt},
    respectively. */

void hr3_L_inf_normalize_point(hr3_point_t *p);
void hr3_L_inf_normalize_plane(hr3_plane_t *P);
void hr3_L_inf_normalize_line(hr3_line_t *L);
  /* Scales the homogeneous coordinates of the given object by a positive
    factor so that the maximum absolute value among them is 1.
    If all coordinates are zero, they are turned into NaNs. */
  
hr3_point_t hr3_point_throw(void);
  /* Returns a point with homogeneous coordinates uniformly
    distributed over the unit sphere of {\RR^4}. */
   
hr3_plane_t hr3_plane_throw(void);
  /* Returns a plane with homogeneous coefficienys uniformly
    distributed over the unit sphere of {\RR^4}. */

/* PRINTOUT 

  The procedures below print the string {pre} (if not NULL) before the object, and 
  and {suf} (if not NULL) after it.  The coordinates/coefficients are printed
  with the format {fmt} (or with "%24.16e" if {fmt} is NULL. */

void hr3_point_print(FILE *wr, char *pre, hr3_point_t *p, char *fmt, char *suf);
  /* Prints {p} to {wr} as "[ {w} {x} {y} {z} ]", prefixed by {pre}
    and followed by {suf}, formatting  each coordinated with {fmt}. */
    
void hr3_plane_print(FILE *wr, char *pre, hr3_plane_t *P, char *fmt, char *suf);
  /* Prints {P} to {wr} as "< {W} {X} {Y} {Z} >", prefixed by {pre} 
    and followed by {suf}, formatting each coefficient with {fmt}. */
    
void hr3_line_print(FILE *wr, char *pre, hr3_line_t *L, char *fmt, char *suf);
  /* Prints {L} to {wr} as "[[ {wx} {wy} {xy} {wz} {xz} {yz} ]]", 
    prefixed by {pre} and followed by {suf},
    formatting each Grassmann coordinate with {fmt}. */

#endif
