/* Oriented projective geometry in three dimensions. */
/* Last edited on 2023-10-09 21:03:08 by stolfi */ 

#ifndef hr3_H
#define hr3_H

/* Based on HR3.i3, created 1993-04-18 by Marcos C. Carrard. */
/* Based on H3.pas by J. Stolfi. */
   
#define _GNU_SOURCE

#include <r3.h>
#include <r3x3.h>
#include <r4.h>
#include <r4x4.h>
#include <r6.h>

#include <sign.h>

typedef struct hr3_point_t { r4_t c; } hr3_point_t; /* c[0..3] are coords [w,x,y,z]. */
typedef struct hr3_plane_t { r4_t f; } hr3_plane_t; /* f[0..3] are coeffs <W,X,Y,Z>. */
typedef struct hr3_line_t  { r6_t k; } hr3_line_t;  /* k[0..5] are Plücker coords [wx,wy,xy,wz,xz,yz] */
  
hr3_point_t hr3_from_r3(r3_t *c);
  /* Point on ``hither'' half of space (i.e, with positive weight)
    whose Cartesian coordinates are {c}. */
    
r3_t r3_from_hr3(hr3_point_t *p);
  /* Cartesian coordinates of point {p} (which must be finite). */

double hr3_pt_pt_diff(hr3_point_t *p, hr3_point_t *q);
  /* Distance between {p} and {q} in the spherical model; that is,
     angle between the vectors {p.c} and {q.c} in {R^4}, in radians. */

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
    
hr3_point_t hr3_point_at_infinity(r3_t *dir);
  /* The point at infinity whose direction, as seen from any hither
    point, is {*dir}. */

double hr3_dist(hr3_point_t *a, hr3_point_t *b);
  /* Distance between {a} and {b}, which must lie in the front half-plane. */
    
double hr3_dist_sqr(hr3_point_t *a, hr3_point_t *b);
  /* Distance squared between {a} and {b}, which must lie in the front half-plane. */
    
r3_t hr3_line_dir(hr3_line_t *L);
  /* The direction along line {L}.  Assumes {L} is not at infinity. */
    
r3_t hr3_plane_normal(hr3_plane_t *P);
  /* The normal direction of plane {P}, on the hither side, pointing into 
    {P}'s positive halfspace.  Assumes {P} is not at infinity. */
    
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

/* PROJECTIVE MAPS */

typedef struct hr3_pmap_t { r4x4_t dir; r4x4_t inv; } hr3_pmap_t;  
  /* A projective map. Field {dir} is the map's matrix, {inv} is its inverse. */

bool_t hr3_pmap_is_identity(hr3_pmap_t *M);
  /* TRUE iff {M} is the identity map (apart from homogeneous scaling). */

hr3_point_t hr3_pmap_point(hr3_point_t *p, hr3_pmap_t *M);
  /* Applies projective map {M} to point {p}. */

hr3_point_t hr3_pmap_inv_point(hr3_point_t *p, hr3_pmap_t *M);
  /* Applies the inverse of projective map {M} to point {p}. */

hr3_plane_t hr3_pmap_plane(hr3_plane_t *P, hr3_pmap_t *M);
  /* Applies projective map {M} to plane {P} */

hr3_plane_t hr3_pmap_inv_plane(hr3_plane_t *P, hr3_pmap_t *M);
  /* Applies the inverse of projective map {M} to plane {P} */

r3_t hr3_pmap_r3_point(r3_t *p, hr3_pmap_t *M);
  /* Applies projective map {M} to the point with Cartesian coordinates {p}. 
    If the result is at infinity, returns {(NAN,NAN,NAN)}. */

r3_t hr3_pmap_inv_r3_point(r3_t *p, hr3_pmap_t *M);
  /* Applies the inverse of projective map {M} to point with Cartesian {p}. 
    If the result is at infinity, returns {(NAN,NAN,NAN)}. */

hr3_pmap_t hr3_pmap_compose(hr3_pmap_t *M, hr3_pmap_t *N);
  /* Returns the composition of {M} and {N}, applied in that order */
    
hr3_pmap_t hr3_pmap_inv(hr3_pmap_t *M);
  /* Returns the inverse of map {M}. */
  
hr3_pmap_t hr3_pmap_inv_compose(hr3_pmap_t *M, hr3_pmap_t *N);
  /* Returns the composition of the inverse of {M} and {N}, applied in that order. */

hr3_pmap_t hr3_pmap_translation(r3_t *v);
  /* Returns the projective map {M} that performs an Euclidean
    translation by the Cartesian vector {v}. */
    
hr3_pmap_t hr3_pmap_scaling(r3_t *scale);
  /* Returns a projective map that performs a scaling of each Cartesian
    coordinate {j} by the factor {scale.c[j]} (that should not be zero).
    The map is actually linear map of {\RR^3} (has {[1 0 0 0]} as the
    first column and the first row). */

hr3_pmap_t hr3_pmap_u_v_rotation(r3_t *u, r3_t *v);
  /* Returns the projective map {M} that performs an Euclidean
    rotation, around some axis through the origin, that takes the
    Cartesian unit vector {u} to the Cartesian unit vector {v} by the
    shortest route.
    
    If {u} is equal to {v}, it returns the identity map. If {u} is
    opposite to {v}, returns a map that performs a 180 degree rotation
    around a random axis through the origin that is orthogonal to
    both. */

hr3_pmap_t hr3_pmap_aff_from_mat_and_disp(r3x3_t *E, r3_t *d);
  /* Returns an affine map {M} that performs the linear map 
    described by the matrix {E} followed by translation by {d}.
    That is, maps the point with Cartesian coordinates {p \in \RR^3} 
    to {p*E + d}. The map will have unit weight (that is, 
    {M.dir.c[0][0] == 1}). */

hr3_pmap_t hr3_pmap_aff_from_four_points(r3_t *o, r3_t *p, r3_t *q, r3_t *r);
  /* Returns a projective map that takes the Cartesian points {(0,0,0)},
    {(1,0,0)}, {(0,1,0)}, and {(0,0,1)} to {o}, {p}, {q}, and {r}, respectively.
    
    The procedure returns a degenerate (non-invertble) map
    if {o,p,q,r} are coplanar.  */

hr3_pmap_t hr3_pmap_from_five_points(hr3_point_t *p, hr3_point_t *q, hr3_point_t *r, hr3_point_t *s, hr3_point_t *u);
  /* Returns a projective map that takes the cardinal points {[1,0,0,0]},
    {[0,1,0,0]}, {[0,0,1,0]}, and {[0,0,0,1]} to {p}, {q}, {r}, and {s}, respectively; and
    also some point {t} of the form {[±1,±1,±1]} to {u}. 
    
    The point {t} is unique, and is the signature of {u} relative to
    the ordered quadruple {p,q,r,s}.
    
    The procedure fails if the set {p,q,r,s,u} contains four coplanar
    points.  */

hr3_pmap_t hr3_pmap_persp(hr3_point_t *obs, hr3_point_t *foc, double rad, hr3_point_t *upp);
  /* Computes a perspective transformation with given viewing parameters:
      {obs} == the scenesys coords of imagesys (0,0,d) (the observer);
      {foc} == the scenesys coords of imagesys (0,0,0) (the image focus);
      {upp}  == the scenesys coords of some point 
               with imagesys X == 0, Y > 0 (up reference).
    
    The resulting projective map takes the observer to the
    {Z}-infinity point {(0,0,+oo)}, and the point {ctr} to the origin.
    Then perspective projection reduces to applying the map and
    discarding the {Z} coordinate.
    
    The virtual camera is tilted around the projection axis {L} (the
    line from {obs} to the origin) in such a way that point {upp} will
    project somewhere onto the positive {Y} axis.
    
    If {rad > 0}, the transformation is then composed with a uniform
    scaling by {1/rad}. Thus the circle of radius {rad} centered at
    {ctr} and normal to the line {obs--ctr} is mapped to the unit
    circle of the XY plane. */
  
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

void hr3_point_print (FILE *f, char *pre, hr3_point_t *a, char *fmt, char *suf);
  /* Prints point {a} to file {f} as "[{w} {x} {y} {z}]". */

void hr3_plane_print (FILE *f, char *pre, hr3_plane_t *P, char *fmt, char *suf);
  /* Prints plane {P} to file {f} as "<{W} {X} {Y} {Z}>". */

void hr3_line_print (FILE *f, char *pre, hr3_line_t *L, char *fmt, char *suf);
  /* Prints line {L} to file {f} as "[{wx} {wy} {xy} {wz} {xz} {yz}]". */

void hr3_pmap_print (FILE *wr, hr3_pmap_t *M, char *pref, char *suff);
  /* Prints {M} on file {wr}, with some default format.  The printout
    starts with the given {pref}, if not {NULL}, and ends with the given {suff},
    if not {NULL}. */

void hr3_pmap_gen_print
  ( FILE *wr, hr3_pmap_t *M,
    char *fmt, char *pref,                /* Overall prefix. */
    char *rpref, char *rsep, char *rsuff, /* Row prefix, matrix separator, and suffix. */
    char *elp, char *esep, char *erp,     /* Delimiters for each row. */
    char *suff                            /* Overall sufffix. */
  );
  /* Prints the projective map {M} to file {wr}, using {fmt} for each
    matrix element.  The printout consists of the string {pref},
    followed by one section for each row of the matrices {M.dir} and {M.inv}, followed by the string {suff}.
    Each section has the string {rpref}, a row of the matrix {M.dir},
    the string {rsep}, a row of matrix {M.inv}, and the string {rsuff}.
    Each row of each matrix is bounded by {elp} and {erp}, and elements are
    separated by {esep}. Defaults are provided for any of these strings which are
    NULL. */

#endif
