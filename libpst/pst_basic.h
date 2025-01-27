#ifndef pst_basic_H
#define pst_basic_H

/* pst_basic.h -- basic data types for gauge-based photostereo. */
/* Last edited on 2025-01-22 18:52:25 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <values.h>

#include <float_image.h>
#include <vec.h>
#include <r2.h>
#include <r3.h>
#include <frgb.h>
#include <argparser.h>

#define INF INFINITY
  /* Plus infinity. */

vec_typedef(image_vec_t,image_vec,float_image_t *);
  /* Defines the type {image_vec_t} as a vector of {float_image_t*} elems. */

typedef r3_t pst_normal_func_t (r2_t *p);
  /* A procedure that computes the normal direction {nrm} at a visible point
    {P} of a some surface, given the projection {p} of that point in some
    plane. 
    
    Both {p} and {nrm} are given in some orthogonal {U,V,W} coordinate
    system such that the projection of point {(u,v,w)} has coordinates
    {(u,v)} (i.e., such that the {W} axis is parallel to the direction
    of projection). The returned normal should have a non-negative {W}
    component.
    
    The procedure should return {(NAN,NAN,NAN)} if 
    the normal direction is not defined at the point
    {P} (e.g. if {P} is at infinity).  Otherwise it should 
    return a valid unit-length vector. */ 
 
typedef frgb_t pst_albedo_func_t (r2_t *p);
  /* A procedure that computes the albedo (intrinsic color) {alb} at a
    visible point {P} of a some surface, given the projection {p} of
    that point in some plane. The components of {alb} must be finite and
    usually (but not necessarily) between 0 and 1.  
    
    The procedure should return the invalid color
    {(NAN,NAN,NAN)} if the albedo is not defined at {p}.
    Otherwise it should return an {frgb_t} color 
    with finite components, usually (but not necessarily)
    between 0 and 1. */ 

#endif
