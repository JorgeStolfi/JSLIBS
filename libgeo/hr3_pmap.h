/* Oriented projective maps in three dimensions. */
/* Last edited on 2025-01-05 00:25:59 by stolfi */ 

#ifndef hr3_pmap_H
#define hr3_pmap_H
   
#include <stdint.h>

#include <sign.h>
#include <r3.h>
#include <r3x3.h>
#include <r4.h>
#include <r4x4.h>
#include <r6.h>
#include <hr3.h>


typedef struct hr3_pmap_t { r4x4_t dir; r4x4_t inv; } hr3_pmap_t;  
  /* A projective map. Field {dir} is the map's matrix, {inv} is its inverse. */

bool_t hr3_pmap_is_identity(hr3_pmap_t *M);
  /* TRUE iff {M} is the identity map (apart from homogeneous scaling). */

hr3_point_t hr3_pmap_point(hr3_point_t *p, hr3_pmap_t *M);
hr3_point_t hr3_pmap_inv_point(hr3_point_t *p, hr3_pmap_t *M);
  /* Applies projective map {M} or its inverse, respectively, to point {p}. */

r3_t hr3_pmap_r3_point(r3_t *p, hr3_pmap_t *M);
r3_t hr3_pmap_inv_r3_point(r3_t *p, hr3_pmap_t *M);
  /* Maps the Cartesian point {p} by the projective map {M} or its inverse,
    respectively, and returns the result converted to Cartesian coordinates. 
    If the result is at infinity, returns {(NAN,NAN,NAN)}. */

hr3_plane_t hr3_pmap_plane(hr3_plane_t *P, hr3_pmap_t *M);
hr3_plane_t hr3_pmap_inv_plane(hr3_plane_t *P, hr3_pmap_t *M);
  /* Applies projective map {M} or its inverse, respectively, to plane {P} */

hr3_pmap_t hr3_pmap_compose(hr3_pmap_t *M, hr3_pmap_t *N);
  /* Returns the composition of {M} and {N}, applied in that order */
    
hr3_pmap_t hr3_pmap_inv(hr3_pmap_t *M);
  /* Returns the inverse of map {M}. */
  
hr3_pmap_t hr3_pmap_inv_compose(hr3_pmap_t *M, hr3_pmap_t *N);
  /* Returns the composition of the inverse of {M} and {N}, applied in that order. */
    
bool_t hr3_pmap_is_affine(hr3_pmap_t *M);
  /* Returns true iff {M} is an affine map; that is, the first column of
    its direct matrix is {w,0,0,0} for some positive {w}. Note that a
    proper projective map is affine if and only if its inverse is
    affine. */

double hr3_pmap_diff_sqr(hr3_pmap_t *M, hr3_pmap_t *N);
  /* Computes the sum of squares differences between corresponding
    elements of {M.dir} and {N.dir} and of {M.inv} and {N.inv},
    after implicitly scaling both so that the sum of squared elements is 1. */

double hr3_pmap_mismatch_sqr(hr3_pmap_t *M, uint32_t np, r3_t p1[], r3_t p2[]);
  /* Computes the mean squared distance between the
     mapped points {p1[0..np-1]} mapped by {M.dir} and the 
     points {p2[0..np-1]} mapped by {M.inv}. */

double hr3_pmap_deform_sqr(r3_t ph[], hr3_pmap_t *M);
  /* Measures the amount of deformation produced by the projective map
    {M} on the cuboid {Q} whose corners {Q[0..1,0..1,0..1]} are stored 
    in {ph[0..7]}, linearized. For meaningful results, the points should
    be in general position (no four of them collinear).
    
    The procedure compares each of the 12 sides and each of the 4 main
    diagonals of {Q} with those of the cuboid {Q'} that is {Q} mapped by
    {M}. For each of these line segments, it computes the original
    length {d} and the length {d'} of the mapped segment, and the log of
    the ratio {d'/d}. The result is the variance of these logs.
    
    Thus the result will be zero if and an only if every distance {d'}
    is the corresponding distance {d} times a common factor {s}; that
    is, if {Q} is similar to {Q'}, meaning that {M} is an similarity
    -- a translation combined with a rotation or mirroring and a uniform
    change of scale (by the common factor {s}). */

double hr3_pmap_aff_discr_sqr(hr3_pmap_t *M, hr3_pmap_t *N);
  /* Assumes that {M} and {N} are affine maps.  Returns the total
    squared mismatch between them, defined as 
    
      { INTEGRAL { |(A - B)(u)|^2 du: u \in \RS^2 } }
    
    */

/* SPECIAL PROJECTIVE MAPS */

hr3_pmap_t hr3_pmap_identity(void);
  /* Returns the identity projective map, defined by the iedntity 4x4 matrix. */

hr3_pmap_t hr3_pmap_translation(r3_t *v);
  /* Returns the projective map {M} that performs an Euclidean
    translation by the Cartesian vector {v}. */
    
hr3_pmap_t hr3_pmap_scaling(r3_t *scale);
  /* Returns a projective map that performs a scaling of each Cartesian
    coordinate {j} by the factor {scale.c[j]} (that should not be zero).
    The map is actually linear map of {\RR^3} (has {[1 0 0 0]} as the
    first column and the first row). */

hr3_pmap_t hr3_pmap_u_to_v_rotation(r3_t *u, r3_t *v);
  /* Returns the projective map {M} that performs an Euclidean
    rotation, around some axis through the origin, that takes the
    Cartesian unit vector {u} to the Cartesian unit vector {v} by the
    shortest route.  See {r3x3_u_to_v_rotation} for details. */

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

/* PRINTOUT */

void hr3_pmap_print(FILE *wr, hr3_pmap_t *M, char *pref, char *suff);
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
