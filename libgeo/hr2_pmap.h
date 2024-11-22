/* Oriented projective maps in two dimensions. */
/* Last edited on 2024-11-21 21:17:19 by stolfi */ 
   
#ifndef hr2_pmap_H
#define hr2_pmap_H

#define _GNU_SOURCE
#include <stdint.h>

#include <sign.h>
#include <r2.h>
#include <r2x2.h>
#include <r3.h>
#include <r3x3.h>
#include <hr2.h>


typedef struct hr2_pmap_t { r3x3_t dir; r3x3_t inv; } hr2_pmap_t;
  /* A projective map. Field {dir} is the map's matrix, {inv} is its inverse. */
 
bool_t hr2_pmap_is_valid(hr2_pmap_t *M, double tol);
  /* Returns true iff {M} is a projective map. Namely, if the product of
    {M.inv} and {M.dir} is the identity, apart from roundoff errror and
    homogeneous scaling of the whole matrices, and the absolute values of 
    the determinants of both {M.dir} and {M.inv} are at least {tol}. */

hr2_point_t hr2_pmap_point(hr2_point_t *p, hr2_pmap_t *M);
hr2_point_t hr2_pmap_inv_point(hr2_point_t *p, hr2_pmap_t *M);
  /* Applies projective map {M} or its inverse, respectively, to point {p}. */

r2_t hr2_pmap_r2_point(r2_t *p, hr2_pmap_t *M);
r2_t hr2_pmap_inv_r2_point(r2_t *p, hr2_pmap_t *M);
  /* Maps the Cartesian point {p} by the projective map {M} or its inverse,
    respectively, and returns the result converted to Cartesian coordinates. 
    If the result is at infinity, returns {(NAN,NAN)}. */

hr2_line_t hr2_pmap_line(hr2_line_t *L, hr2_pmap_t *M);
hr2_line_t hr2_pmap_inv_line(hr2_line_t *L, hr2_pmap_t *M);
  /* Applies projective map {M} or its inverse, respectively, to line {L}. */

hr2_pmap_t hr2_pmap_compose(hr2_pmap_t *M, hr2_pmap_t *N);
  /* Returns the composition of {M} and {N}, applied in that order. */

hr2_pmap_t hr2_pmap_inv(hr2_pmap_t *M);
  /* Returns the inverse of map {M}. */

hr2_pmap_t hr2_pmap_inv_compose(hr2_pmap_t *M, hr2_pmap_t *N);
  /* Returns the composition of the inverse of {M} and {N}, applied in that order. */

double hr2_pmap_diff_sqr(hr2_pmap_t *M, hr2_pmap_t *N);
  /* Computes the sum of squares differences between corresponding
    elements of {M.dir} and {N.dir} and of {M.inv} and {N.inv},
    after implicitly scaling both so that the sum of squared elements is 1. */

double hr2_pmap_mismatch_sqr(hr2_pmap_t *M, uint32_t np, r2_t p1[], r2_t p2[], double w[]);
  /* Computes the weighted mean squared distance between the
     points {p1[0..np-1]} and the points {p2[0..np-1]} under the map {M}.
     Specifically, it is half of the weighted average of {|M.dir(p1[i])-p2[i]|^2} and
     {|p1[i]-M.inv(p2[i])|^2}, over all {i}, each pair taken with weight {w[i]}.
     The weights must be non-negative.  If {w} is {NULL}, assumes
     all equal weights.  If {np} is zero or the weights are all zero,
     returns 0.0. */

double hr2_pmap_max_mismatch(hr2_pmap_t *M, uint32_t np, r2_t p1[], r2_t p2[]);
  /* Compute max distance between {p1[k]} mapped by {M} and {p2[k]}
    and between {p1[k]} and {p2[k]} mapped by {M^{-1}}, for {k} in {0..np-1}. */

void hr2_pmap_show_point_mismatch(hr2_pmap_t *M, uint32_t np, r2_t p1[], r2_t p2[], double w[]);
  /* Prints to {stderr} the mismatch between {M.dir(p1[k])} and {p2[k]}, and
    between {p1[k]} and {M.inv(p2[k])}, for {k} in {0..np-1}. */

double hr2_pmap_deform_sqr(r2_t ph[], hr2_pmap_t *M);
  /* Measures the amount of deformation produced by the projective map
    {M} on the set of four points {Q} whose corners are {ph[0..3]}, in
    any order.  For meaningful results, the points should be
    in general position (no three collinear).
    
    The procedure compares each of the 6 point-to-point distances 
    of {Q} with those of the set {Q'} that is {Q} mapped by {M}.
    For each of these line segments, it computes the original length {d} 
    and the length {d'} of the mapped segment, and the log of the 
    ratio {d'/d}. The result is is the variance of these logs.
    
    Thus the result will be zero if and an only if every distance {d'}
    is the corresponding distance {d} times a common factor {s}; that
    is, if {Q} is similar to {Q'}, meaning that {M} is an similarity
    -- a translation combined with a rotation or mirroring and a uniform
    change of scale (by the common factor {s}). */

/* SPECIAL PROJECTIVE MAPS */

hr2_pmap_t hr2_pmap_identity(void);
  /* Returns the identity projective map, defined by the iedntity 3x3 matrix. */

hr2_pmap_t hr2_pmap_mirror(sign_t sgnx, sign_t sgny);
  /* Returns a projective map that preserves or negates the {x}
    Cartesian coordinate, depending on {sgnx} being {-1} or {+1};
    and similarly for {sgny} and the {y} coordinate.  The signs must not be zero. */

hr2_pmap_t hr2_pmap_xy_swap(void);
  /* Returns a projective map that swaps the {x} and {y} Cartesian coordinates. */

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

hr2_pmap_t hr2_pmap_from_four_points(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r, hr2_point_t *u);
  /* Returns a projective map that takes the cardinal points {[1,0,0]},
    {[0,1,0]}, and {[0,0,1]} to {p}, {q}, and {r}, respectively; and
    also some point {s} of the form {[±1,±1,±1]} to {u}. 
    
    The point {s} is unique, and is the signature of {u} relative to
    the ordered triple {p,q,r}.
    
    The procedure fails if the set {p,q,r,u} contains three collinear
    points. */

hr2_pmap_t hr2_pmap_from_four_r2_points(r2_t *p, r2_t *q, r2_t *r, r2_t *u);
  /* Returns a projective map that takes the hither points with Cartesian coordinates {(0,0)},
    {(1,0)}, and {(0,1)} to points with Cartesian coordinates {p}, {q},
    and {r}, respectively; and also some point {s} to the point {u}.
    
    The point {s} is one of four points of the form {(+1,+1)}, {(-1,+1)},
    {(+1,-1)}, or {(1/3,1/3)}, depending on the position of {u} 
    relative to the ordered triple {p,q,r}.
    
    The procedure fails if the set {p,q,r,u} contains three collinear
    points. */

hr2_pmap_t hr2_pmap_from_parms
  ( r2_t *persp,
    double shear,
    double skew,
    double scale,
    double ang,
    r2_t *disp
  );
  /* Returns a projective map with the given parameters.  
  
    The map is a composition of a perspective map {P(persp)} and three
    affine maps: a shear and non-uniform scaling map {S(shear,skew)}, a
    rotation and uniform scaling map {R(scale,ang)}, and a translation
    map {T(disp)}.
    
    The map {P(persp)} will map each Cartesian point {(X,Y)} to
    {(X,Y)/(1 + a*X +b*Y)} where {(a,b)} is the vector {persp}. The is
    map keeps fixed the origin and all points on the line through the
    origin perpendicular to {(a,b)}. Its Jacobian at the origin is the
    {2×2} identity - that is, infinitesimally small figures at the
    origin are not changed.  If {persp} is {NULL}, assumes {P((0,0))} 
    (the identity).
    
    The map {S(shear,skew)} keeps the origin fixed, maps the {\RR^2}
    points {(1,0)} and {(0,1)} respectively to {(1,shear)*rsk} and
    {(shear,1)/rsk}, where {rsk = sqrt(skew)}. In particular,
    {S(0,skew)} scales {X} by {rsk} and {Y} by {1/rsk}. It requires
    {|shear| < 1} and {skew > 0}.
    
    The map {R(scale,ang)} keeps the origin fixed, scales everything by
    {scale}, and rotates the plane by {ang} about the origin. It
    requires {scale > 0}.
    
    The map {T(disp)} displaces evry point by the vector {disp}.
    If {disp} is {NULL}, assumes {T((0,0))} (the identity).
    
    The affine map will have {M.dir[0,0] = M.inv[0,0] = 1.0}. */


void hr2_pmap_throw(hr2_pmap_t *M);
  /* Fills {M} with a random projective map with definitely nonzero 
    determinant. */

hr2_pmap_t hr2_pmap_r2_from_sign_class(uint32_t class);
  /* Returns a projective map that takes the hither points with Cartesian 
    coordinates {(0,0)}, {(1,0)}, and {(0,1)} to points with the same
    Cartesian coordinates.  
    
    The weight of the images depend on {class},  an integer in {0..3]}.
    The point {[1,0,0]} is always mapped to itself.  The image of {[1,1,0]}
    will be {[1,1,0]} or {[-1,-1,0]} depending on bit {class&1} is 0 or 1,
    respectively.  Likewise, the image of {[1,0,1]} will be {[1,0,1]} or {[-1,0,-1]}
    depending on bit {class&2}.
    
    The image of the barycenter {[3,1,1]} of the three source points will be
    a point {u} whose Cartesian coordinates depend on the {class}, an integer in {0..3}.
    Namely,
    
      {class = 0} --> {u = (1/3,1/3)}
      {class = 1} --> {u = (+1,-1)}
      {class = 2} --> {u = (-1,+1)}
      {class = 3} --> {u = (+1,+1)}
      
    In fact, if {class = 0}, the map is the identity.  In the other cases,
    the map takes some hither points to the yonder side, and vice-versa. */
  
sign_t hr2_pmap_sign(hr2_pmap_t *M);
  /* Returns {+1}, 0, or {-1} depending on whether the determinant of {M}
    is positive, zero, or negative. 
    
    Note that the if the determinant is actually zero or close to zero
    the result may be any of the three, because of floating-point
    roundoff errors. */
  
/* MAP TYPE CHECKING

  The following procedures check whether the map {M} is of the type
  specified by {type} and sign parameters.  All these procedures
  allow from homogeneous scaling of the matrices, and noise of
  magnitude {M[0][0]*tol} in all elements. */
    
/* SPECIAL MAP TYPES */

typedef enum 
  { hr2_pmap_type_IDENTITY,
    hr2_pmap_type_TRANSLATION,
    hr2_pmap_type_CONGRUENCE,
    hr2_pmap_type_SIMILARITY,
    hr2_pmap_type_AFFINE,
    hr2_pmap_type_GENERIC,
    hr2_pmap_type_NONE
  } hr2_pmap_type_t;
  /* Identifies a type of projective map.  
  
    A congruence is a combination
    of translation and rotation or reflection.  A similarity is the same 
    plus a uniform scaling. An affine map is a linear map of {\RR^2}
    combined with a translation.  The value {hr2_pmap_type_GENERIC} means
    a general (unrestricted) projecive map. T
    
    he value {hr2_pmap_type_NONE} (which must be the last in the enumeration)
    is a null value that can be used to say for "unspecified" or "invalid". */

#define hr2_pmap_type_FIRST hr2_pmap_type_IDENTITY
  /* The first value in the enum {hr2_pmap_type_t}. */

#define hr2_pmap_type_LAST hr2_pmap_type_GENERIC
  /* The last value in the enum {hr2_pmap_type_t} other than {hr2_pmap_type_NONE}. */

bool_t hr2_pmap_is_identity(hr2_pmap_t *M, double tol);
  /* TRUE iff {M} is the identity map, apart from noise of 
    relative magnitude {tol}. */

bool_t hr2_pmap_is_translation(hr2_pmap_t *M, double tol);
  /* TRUE iff {M} is a translation, apart from noise of 
    relative magnitude {tol}. */
   
bool_t hr2_pmap_is_congruence(hr2_pmap_t *M, double tol);
  /* TRUE iff {M} is a congruence, apart from noise of 
    relative magnitude {tol}). */

bool_t hr2_pmap_is_similarity(hr2_pmap_t *M, double scale, double tol);
  /* TRUE iff {M} is a similarity (apart from homogeneous scaling and
    noise of magnitude {tol}).  If {scale} is not {NAN}, requires the
    scaling factor to be {scale} (which must be positive). */

bool_t hr2_pmap_is_affine(hr2_pmap_t *M, double tol);
  /* Returns true iff {M} is a POSITIVE affine map (apart from homogeneous
    scaling and noise of magnitude {tol}). In particular, returns
    {FALSE} if the determinant of {M.dir} or {M.inv} is negative
    or less than {tol}.  */

bool_t hr2_pmap_type_from_string(char *tname, hr2_pmap_type_t *type_P);
  /* Converts the string {tname} to a {hr2_pmap_type_t} value.  If it succeeds,
    it stores the type into {*type_P} and returns {TRUE}. 
    If the string is not a valid type, leaves {*typ_P} unchanged and
    returns {FALSE}.
    
    Specifically, maps "IDENTITY" and "identity" to {hr2_pmap_type_IDENTITY},
    "TRANSLATION" and "translation" to {hr2_pmap_type_TRANSLATION},
    and so on. */

char *hr2_pmap_type_to_string(hr2_pmap_type_t type);
  /* Converts the {type} to a string, all in upper case. That is, maps
    {hr2_pmap_type_IDENTITY} to "IDENTITY", {hr2_pmap_type_TRANSLATION}
    to "TRANSLATION", and so on. The returned string should be treated
    as read-only and should not be given to {free}. */
    
void hr2_pmap_invert_sign(hr2_pmap_t *M);
  /* Inverts the handedness of {M} by swapping rows 
    1 and 2 of the direct matrix, and columns 1 and 2 of the 
    inverse matrix. */

void hr2_pmap_set_sign(hr2_pmap_t *M, sign_t sgn);
  /* Makes sure that {M} has the handedness {sgn} (either {-1} or {+1}), by
    applying {hr2_pmap_invert_sign}.  Fails if the determinant
    of {M} is exactly zero. */

bool_t hr2_pmap_is_type(hr2_pmap_t *M, hr2_pmap_type_t type, sign_t sgn, double tol);
  /* If {sgn} is {+1}, returns true iff the map {*M} is of the given
    type (apart from homogeneous scaling) and positive determinant, with
    relative tolerance {tol}.
    
    If {sgn} is {-1}, {M} should have negative determinant. In that
    case, returns true iff {M} is the negative version of a map of the
    given type.
    
    If {sgn} is zero, returns true if either of the above conditions is true.
    
    Returns false if the appropriate condition is not satisfied. In
    particular, if the relative magnitude of the determinant of either
    {M.dir} or {M.inv} is less than {tol}.*/

void hr2_pmap_print (FILE *wr, hr2_pmap_t *M, char *pref, char *suff);
  /* Prints {M} on file {wr}, with some default format.  The printout
    starts with the given {pref}, if not {NULL}, and ends with the given {suff},
    if not {NULL}. */

void hr2_pmap_gen_print
  ( FILE *wr, hr2_pmap_t *M,
    char *fmt,                            /* Format for elements. */
    char *pref,                           /* Overall prefix. */
    char *rpref, char *rsep, char *rsuff, /* Row prefix, matrix separator, and row suffix. */
    char *elp, char *esep, char *erp,     /* Delimiters for matrix in each row. */
    char *suff                            /* Overall suffix. */
  );
  /* Prints the projective map {M} to file {wr}, using {fmt} for each
    matrix element.  
    
    The printout consists of the string {pref}, followed by three rows
    with the matrices {M.dir} and {M.inv}, followed by the string
    {suff}. Each row of the matrices has the string {rpref}, a row of
    the matrix {M.dir}, the string {rsep}, a row of matrix {M.inv}, and
    the string {rsuff}. Each row of each matrix is bounded by {elp} and
    {erp}, and elements are separated by {esep}. Defaults are provided
    for any of these strings which are NULL. */

#endif
