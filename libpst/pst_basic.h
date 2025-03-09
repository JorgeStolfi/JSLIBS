#ifndef pst_basic_H
#define pst_basic_H

/* pst_basic.h -- basic data types for gauge-based photostereo. */
/* Last edited on 2025-03-01 19:53:05 by stolfi */

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

void pst_perturb_normal(r3_t *nrm, double amt);
  /* Randomly perturbs the unit vector {nrm} by the amount {amt},
    preserving its unit length. In particular, {amt=0} means no
    perturbation, while {amt=1} means that the result is a uniformly
    distributed on the unit sphere and independent of the given {nrm}.
    For small {amt}, the perturbation is (to first order) a tangential
    vector perpendicular to {nrm}, with zero mean and root mean square
    length {amt}.  */

r2_t pst_slope_from_normal(r3_t *nrm, double maxSlope);
  /* Converts a normal vector {nrm} to a gradient vector (that is, the
    slopes {dZ/dX} and {dZ/dY}). 
    
    The Z coordinate or the normal must be non-negative. The length of
    the normal vector needs not be 1. If necessary, the computed
    gradient is scaled so that its Euclidean norm does not exceed
    {maxSlope}. The result is undefined if {nrm} is the zero
    vector. */

r3_t pst_normal_from_slope(r2_t *grd);
  /* Converts a gradient vector {grd} (that is, the slopes {dZ/dX} and
    {dZ/dY}) to an outwards-pointing normal vector, with unit
    Euclidean length. */
    
void pst_ensure_pixel_consistency(int32_t NC, int32_t wch, bool_t bad, float v[]);
  /* Assumes that {v[0..NC-1]} are the samples of a pixel
    of some image, typicall a height, slope, or normal map;
    and that {v[wch]}, if {wch} is in {0..NC-1}, is the 
    reliability weight for that pixel, which must be finite an 
    non-negative.
    
    If {bad} is true, or {v[wch]} exists but is zero, or any of the
    samples {v[0..NC-1]} is not finite, sets {v[wch]} (if it exists) to
    0 and all other samples to {NAN}. Otherwise does not change
    {v[0..NC-1]}. */

#endif
