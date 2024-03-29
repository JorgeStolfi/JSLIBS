/* Oriented projective geometry in two dimensions. */
/* Last edited on 2023-10-21 07:03:16 by stolfi */ 
   
#ifndef hr2_H
#define hr2_H

#define _GNU_SOURCE

#include <r2.h>
#include <stdint.h>
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
   
hr2_point_t hr2_point_throw(void);
  /* Returns a point with homogeneous coordinates uniformly
    distributed over the unit sphere of {\RR^3}. */
   
hr2_line_t hr2_line_throw(void);
  /* Returns a line with homogeneous coefficienys uniformly
    distributed over the unit sphere of {\RR^3}. */
 
/* PROJECTIVE MAPS */

typedef struct hr2_pmap_t { r3x3_t dir; r3x3_t inv; } hr2_pmap_t;
  /* A projective map. Field {dir} is the map's matrix, {inv} is its inverse. */

bool_t hr2_pmap_is_identity(hr2_pmap_t *M);
  /* TRUE iff {M} is the identity map (apart from homogeneous scaling). */

hr2_point_t hr2_pmap_point(hr2_point_t *p, hr2_pmap_t *M);
hr2_point_t hr2_pmap_inv_point(hr2_point_t *p, hr2_pmap_t *M);
  /* Applies projective map {M} or its inverse, respectively, to point {p}. */

r2_t hr2_pmap_r2_point(r2_t *p, hr2_pmap_t *M);
r2_t hr2_pmap_inv_r2_point(r2_t *p, hr2_pmap_t *M);
  /* Maps the Cartesian point {p} by the projective map {M} or its inverse,
    respectively, and returns the result converted to Cartesian coordinates. 
    If the result is at  infinity, returns {(NAN,NAN)}. */

hr2_line_t hr2_pmap_line(hr2_line_t *L, hr2_pmap_t *M);
hr2_line_t hr2_pmap_inv_line(hr2_line_t *L, hr2_pmap_t *M);
  /* Applies projective map {M} or its inverse, respectively, to line {L}. */

hr2_pmap_t hr2_pmap_compose(hr2_pmap_t *M, hr2_pmap_t *N);
  /* Returns the composition of {M} and {N}, applied in that order. */

hr2_pmap_t hr2_pmap_inv(hr2_pmap_t *M);
  /* Returns the inverse of map {M}. */

hr2_pmap_t hr2_pmap_inv_compose(hr2_pmap_t *M, hr2_pmap_t *N);
  /* Returns the composition of the inverse of {M} and {N}, applied in that order. */
    
bool_t hr2_pmap_is_affine(hr2_pmap_t *M);
  /* Returns true iff {M} is an affine map; that is, the first column of
    its direct matrix is {w,0,0} for some positive {w}. Note that a
    proper projective map is affine if and only if its inverse is
    affine. */

double hr2_pmap_diff_sqr(hr2_pmap_t *M, hr2_pmap_t *N);
  /* Computes the sum of squares differences between corresponding
    elements of {M.dir} and {N.dir} and of {M.inv} and {N.inv},
    after implicitly scaling both so that the sum of squared elements is 1. */

double hr2_pmap_mismatch_sqr(hr2_pmap_t *M, int32_t np, r2_t p1[], r2_t p2[]);
  /* Computes the mean squared distance between the
     mapped points {p1[0..np-1]} mapped by {M.dir} and the 
     points {p2[0..np-1]} mapped by {M.inv}. */

double hr2_pmap_deform_sqr(r2_t h[], hr2_pmap_t *M);
  /* Measures the amount of deformation produced by the projective map
    {M} on the quadrilateral {Q} whose corners are {h[0..3]}. The
    procedure compares each side and diagonal of {Q} with those of the
    quadrilateral {Q'} that is {Q} mapped by {M}. For each of these line
    segments, it computes the original length {d} and the length {d'} of
    the mapped segment, and the log of the ratio {d'/d}. The result is
    is the variance of these logs.
    
    Thus the result will be zero if and an only if every distance {d'}
    is the corresponding distance {d} times a common factor {s}; that
    is, if {Q} is congruent to {Q'}, meaning that {M} is an similarity
    -- a translation combined with a rotation or mirroring and a uniform
    change of scale (by the common factor {s}). */

double hr2_pmap_aff_discr_sqr(hr2_pmap_t *M, hr2_pmap_t *N);
  /* Assumes that {M} and {N} are affine maps.  Returns the total
    squared mismatch between them, defined as 
    
      { INTEGRAL { |(A - B)(u(t))|^2 : t \in [0 _ 1] } }
    
    where {u(t)} is the unit vector {(cos(2\pi t),sin(2\pi t))}. */

/* SPECIAL PROJECTIVE MAPS */

hr2_pmap_t hr2_pmap_identity(void);
  /* Returns the identity projective map, defined by the iedntity 3x3 matrix. */

hr2_pmap_t hr2_pmap_translation(r2_t *vec);
  /* Returns a projective map that performs a translation by the Cartesian vector {vec}.
    The map is actually affine (has {[1 0 0]} as the first column). */

hr2_pmap_t hr2_pmap_rotation(double ang);
  /* Returns a projective map that performs a rotaton by {ang} radians
    about the origin. The map is actually a linear map of {\RR^2} (has
    {[1 0 0]} as the first column and the first row). */

hr2_pmap_t hr2_pmap_scaling(r2_t *scale);
  /* Returns a projective map that performs a scaling of each Cartesian
    coordinate {j} by the factor {scale.c[j]} (that should not be zero).
    The map is actually linear map of {\RR^2} (has {[1 0 0]} as the
    first column and the first row). */

hr2_pmap_t hr2_pmap_rotation_and_scaling(double ang, double scale);
  /* Returns a projective map that performs a rotaton by {ang} radians
    combined with uniform scaling by {scale}, both about the origin.
    The map is actually linear map of {\RR^2} (has {[1 0 0]} as the
    first column and the first row). */

hr2_pmap_t hr2_pmap_congruence_from_point_and_dir(r2_t *p, r2_t *u, bool_t flip);
  /* Returns a projective map that is actually a Cartesian congruence
    taking the origin {[1,0,0]} to the Cartesian point {p}
    the direction of the {X}-axis to the direction vector {u}.
    The length of {u} is ignored, but must not be zero.
    
    A Cartesian congruence (or isometry) is a map of {\RR^2} to {\RR^2}
    that preserves all distances, and therefore also all angles. It is a
    special case of similarity and of an affine map. It preserves the sign
    of homogeneous coordinate 0 (weight). 
    
    If {flip} is false, the map will preserve handedness, i. e. will be
    a combination of a rotation plus a translation. If {flip} is true,
    the map will reverse handedness; it will be a combination
    of rotation, translation, and mirroring about a line. */

hr2_pmap_t hr2_pmap_similarity_from_two_points(r2_t *p, r2_t *q, bool_t flip);
  /* Returns a projective map that is a Cartesian similarity taking the points {[1,0,0]},
    and {[1,1,0]} to {p} and {q}, respectively.  The points must be distinct.
    
    A similarity (or Euclidean transformation) is a map from {\RR^2} to
    {\RR^2} that preserves ratios between distances, and therefore
    preserves all angles. It is a special case of affine map. It
    preserves the sign of homogeneous coordinate 0 (weight).
    
    If {flip} is false, the map will preserve handedness, i. e. will be
    the combination of a uniform scaling by a positive factor, a
    rotation, and a translation. If {flip} is true, the map will reverse
    handedness, i.e. will be the combination of those and
    and a reflection about some line. Note that a
    uniform scaling by a negative factor is equivalent to a rotation by
    180 degrees.  */

hr2_pmap_t hr2_pmap_aff_from_mat_and_disp(r2x2_t *E, r2_t *d);
  /* Returns an affine map {M} that performs the linear map 
    described by the matrix {E} followed by translation by {d}.
    That is, maps the point with Cartesian coordinates {p \in \RR^2} 
    to {p*E + d}. The map will have unit weight (that is, 
    {M.dir.c[0][0] == 1}). */
     
hr2_pmap_t hr2_pmap_aff_from_three_points(r2_t *o, r2_t *p, r2_t *q);
  /* Returns a projective map that takes the points {(0,0)},
    {(0,1)}, and {(1,0)} to {o}, {p}, and {q}, respectively.
    
    The procedure returns a degenerate (non-invertble) map
    if {o,p,q} are collinear.  */

hr2_pmap_t hr2_pmap_from_four_points(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r, hr2_point_t *u);
  /* Returns a projective map that takes the cardinal points {[1,0,0]},
    {[0,1,0]}, and {[0,0,1]} to {p}, {q}, and {r}, respectively; and
    also some point {s} of the form {[�1,�1,�1]} to {u}. 
    
    The point {s} is unique, and is the signature of {u} relative to
    the ordered triple {p,q,r}.
    
    The procedure fails if the set {p,q,r,u} contains three collinear
    points.  */

typedef enum 
  { hr2_pmat_type_TRANSLATION,
    hr2_pmat_type_CONGRUENCE,
    hr2_pmat_type_SIMILARITY,
    hr2_pmat_type_AFFINE,
    hr2_pmat_type_PROJECTIVE
  } hr2_pmat_type_t;
  /* Identifies a type of projective map.  A congruence is a combination
    of translation and rotation or reflection.  A similarity is the same 
    plus a uniform scaling. An affine map is a linear map of {\RR^2}
    combined with a translation. {hr2_pmat_type_PROJECTIVE} means
    a general (unrestricted) projecive map. */

#define hr2_pmat_type_LAST hr2_pmat_type_PROJECTIVE
  /* The last value in the enum {hr2_pmat_type_t}. */

void hr2_pmap_print (FILE *wr, hr2_pmap_t *M, char *pref, char *suff);
  /* Prints {M} on file {wr}, with some default format.  The printout
    starts with the given {pref}, if not {NULL}, and ends with the given {suff},
    if not {NULL}. */

void hr2_pmap_gen_print
  ( FILE *wr, hr2_pmap_t *M,
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

/* SPECIAL PROCEDURES FOR AFFINE MAPS */

#endif
