/* Affine maps of the plane. */
/* Last edited on 2020-10-16 01:46:24 by jstolfi */ 
   
#ifndef r2_aff_map_H
#define r2_aff_map_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <r2.h>
#include <r3.h>
#include <r2x2.h>
#include <sign.h>

typedef struct r2_aff_map_t 
  { r2x2_t mat;
    r2_t disp;
  } r2_aff_map_t;
  /* An affine map from {\RR^2} to {\RR^2}, specifically that maps 
   a point {p} to the point {p*mat + disp}. The map is degenerate 
   (non-invertible) iff {mat} is singular. */

double *r2_aff_map_elem_addr(r2_aff_map_t *A, int32_t s);
  /* Returns the address of element {s} of the affine map {A}.
    The integer {s} must be in {0..5}.  For {s} is in {0..3},
    returns the addresses of the elements of {A.mat}, in row-by roe order.
    For {s} in {4..5}, returns the addresses of the elements of {A.disp}. */
     
void r2_aff_map_apply(r2_t *p, r2_aff_map_t *A, r2_t *q);
  /* Stores into {*q} the image of {p} by the map {A}. */
   
void r2_aff_map_invert(r2_aff_map_t *A, r2_aff_map_t *B);
  /* Stores into {*B} the inverse of the affine map {*A}. 
    The matrix {f.mat} had better be invertible. */
  
void r2_aff_map_compose(r2_aff_map_t *A, r2_aff_map_t *B, r2_aff_map_t *C);
  /* Stores into {*C} the composition of the maps {A} and {B},
    applied in that order.  That is, if {A} maps {p} to {q},
    and {B} maps {q} to {r}, then {C} maps {p} to {r}. */

r2_aff_map_t r2_aff_map_translation(r2_t *vec);
  /* Returns a projective map that performs a translation by the Cartesian vector {vec}. */

r2_aff_map_t r2_aff_map_rot_scale(double ang, double scale);
  /* Returns a projective map that performs a rotaton by {ang} radians
    combined with uniform scaling by {scale}, both about the origin */

void r2_aff_map_disp_sqr
  ( r2_aff_map_t *A, 
    r2_aff_map_t *B, 
    r2_aff_map_t *R,
    double *dabs2P,
    double *drel2P
  );
  /* Returns in {*dabs2P} and {*drel2P} (if not NULL)
    the total squared discrepancy between {*A} and 
    {*B}, respectively absolute and relative to the adjustment size {*R}.  
    Namely 
    
      { dabs2 = SUM { (Ae[s] - Be[s])^2 } }
    
      { drel2 = SUM { ((Ae[s] - Be[s])/Re[s])^2 } }
    
    where {Ae[s],Be[s],Re[s]} are all corresponding elements of
    {*A,*B,*C}. excluding those where {Re[s]} is zero. 
    
    Clients of this module may want to add a small multiple of
    this function to the goal function in order to bias the search
    towards the neighborhood of the initial guess. */

double r2_aff_map_mismatch_sqr( r2_aff_map_t *A, r2_aff_map_t *B);
  /* Returns the total squared mismatch between {*A} and 
    {*B}, defined as 
    
      { INTEGRAL { |(A - B)(u(t))|^2 : t \in [0 _ 1] } }
    
    where {u(t)} is the unit vector {(cos(2\pi t),sin(2\pi t))}. */

r2_aff_map_t r2_aff_map_from_points(r2_t *o, r2_t *p, r2_t *q);
  /* Returns a projective map that takes the points {(0,0)},
    {(0,1)}, and {(1,0)} to {o}, {p}, and {q}, respectively.
    
    The procedure returns a degenerate (non-invertble) map
    if {o,p,q} are collinear.  */
    
r2_aff_map_t r2_aff_map_from_point_pairs(int32_t np, r2_t p1[], r2_t p2[], double w[]);
  /* Returns an affine map {M} that best maps the points {p1[0..np-1]}
    to the points {p2[0..np-1]}, in the sense of minimizing the
    mean squared error with weights {w[0..np-1]}.  
    
    The weights had better have positive sum.  If {w == NULL}, assumes
    equal weights for all points.
    
    If {np} is 3 or less, the map is exact. If {np} is zero, the result is
    the identity map.  If {np==1}, the result is a translation. If {np==2}, the result is an
    Euclidean similarity (translation, rotation and scaling). If {np}
    is 3 or more the result is a general affine map. */

void r2_aff_map_print (FILE *f, r2_aff_map_t *A);
  /* Prints {A} on file {f}, with some default format. */

void r2_aff_map_gen_print
  ( FILE *f, r2_aff_map_t *A,
    char *mfmt, char *dfmt,
    char *olp, char *osep, char *orp, /* Outer delimiters. */
    char *ilp, char *isep, char *irp,  /* Inner delimiters. */
    char *dsep  /* Matrix to displacement separator. */
  );
  /* Prints the affine map {A} to file {f}, using {mfmt} for each
    matrix element and {dfmt} for each displacement element.
    The matrix {A.mat} is bounded by {olp} and {orp}, and its
    rows are separated by {osep}. Each row of {A.mat} as well as the
    vector {A.disp} is bounded by {ilp} and {irp}, and elements are
    separated by {isep}. The matrix and displacement are separated by
    {dsep}. Defaults are provided for any of these strings which are
    NULL. */


#endif
