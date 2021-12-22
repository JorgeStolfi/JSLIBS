#ifndef float_image_align_H
#define float_image_align_H

/* Tools for optimizing a vector of points on the plane. */
/* Last edited on 2021-12-17 15:01:22 by stolfi */ 

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <float_image.h>
#include <float_image_align.h>

void float_image_align_single_scale_enum
  ( int ni,                           /* Number of objects to align. */
    i2_t iscale,                      /* Object scaling exponent along each axis. */  
    float_image_align_mismatch_t *f2, /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],                      /* Max alignment adjustment for each object. */
    double tol,                       /* Desired precision. */
    r2_t p[],                         /* (IN/OUT) Corresponding points in each object. */
    double *f2p                       /* (OUT) Mismatch for the computed alignment vector. */
  );
  /* Adjusts an alignment vector {p[0..ni-1]} for some {ni} objects so
    that they are as similar as possible in the neighborhoords of those
    points, as defined by the mismatch function {f2(ni,p,iscale)}.
    
    On input, {p[0..ni-1]}, must be a guess for the optimum alignment.
    On output, {p[0..ni-1]} will be the best alignment found.
    
    Uses exhaustive enumeration of alignment vectors {q} that differ from {p} by
    balanced dsppacement vectors and lie within the ellipsoid whose center is
    the given guess {p} and whose radius vector is {arad}. 
    
    Considers any coordinate {p[i].c[j]} fixed if {arad[i].c[j]} is 
    zero.  Otherwise, {arad[i].c[j]} must be positive.
    
    The probe points {q} will comprise a compact regular lattice with closest
    distance {tol} that includes the initial guess {p}.
    
    The value of {f2(ni,iscale,p)} is returned on {*f2p}. */

#endif
