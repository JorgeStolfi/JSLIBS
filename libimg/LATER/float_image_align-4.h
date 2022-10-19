#ifndef float_image_align_H
#define float_image_align_H

/* Tools for optimizing a vector of points on the plane. */
/* Last edited on 2017-06-05 16:02:48 by stolfilocal */ 

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <float_image.h>

/* 
  ALIGNMENT VECTOR
  
  A (two-dimensional) /alignment vector/ for {ni} objects is a list of points 
  {p[0..ni-1]} of the plane.  It says that point {p[i]} of each object {i}
  corresponds in some sense to object {p[j]} of each other object {j}. 
  
  For example, the objects could be images, and the alignment {p[0..ni-1]} could be 
  claiming that some neighborhood of point {p[i]} of image {i} looks like
  the similar neighborhood of point {p[j]} of image {j}. */

typedef double float_image_align_mismatch_t(int ni, r2_t p[], i2_t iscale); 
  /* Type of a function that evaluates the mismatch of {ni} objects.
    in some neighborhood of points {p[0..ni-1]}, at a specified
    scale of resolution {iscale}.
    
    The mismatch should be computed as if every object is expanded or
    reduced by the scale factor {1/2^iscale.c[j]} along axis {j}, with
    each alignment point {p[i]} scaled by that same amount.
    
    The size of the neighborhood, measured on the scaled object, should
    be the same. That is, measured on the unscaled object, the
    neighborhod should be centered on {p[i]} and have size proportional
    to {2^iscale.c[j]} along axis {j}.
    
    Nevertheless, the function should ideally take about the same time
    whatever the {iscale}. For positive {iscale}, this usually requires
    that each object be reduced with antialiasing before the
    comparison.
    
    The function had better be C2 near the optimum alignment for
    quadratic optimization (see below). */

void float_image_align_single_scale_enum
  ( int ni,                   /* Number of objects to align. */
    i2_t iscale,              /* Object scaling exponent along each axis. */  
    float_image_align_mismatch_t *f2,  /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],              /* Max alignment adjustment for each object. */
    r2_t step[],              /* Adjustment step for each object. */
    r2_t p[],                 /* (IN/OUT) Corresponding points in each object. */
    double *f2p               /* (OUT) Mismatch for the computed alignment vector. */
  );
  /* Adjusts an alignment vector {p[0..ni-1]} for some {ni} objects so
    that they are as similar as possible in the neighborhoords of those
    points, as defined by the mismatch function {f2(ni,p,iscale)}.
    Uses exhaustive enumeration of all valid alignment vectors.
    
    On input, {p[0..ni-1]} must be a guess for the optimum alignment. On
    output, {p[0..ni-1]} will be an improved alignment, which differs
    from the given one by an integer multiple of {step[i].c[j]} in each
    axis {j=0..1}, not exceeding {arad[i].c[j]} in absolute value.
    
    If {arad[i].c[j]} is positive, {step[i].c[j]} must be positive. If
    {arad[i].c[j]} is zero or less than {step[i].c[j]}, the coordinate
    {p[i].c[j]} is fixed.
    
    The value of {f2(ni,iscale,p)} is returned on {*f2p}. */

void float_image_align_single_scale_quadopt
  ( int ni,                   /* Number of objects to align. */
    i2_t iscale,              /* Object scaling exponent along each axis. */  
    float_image_align_mismatch_t *f2,  /* Function that evaluates the mismatch between the objects. */
    r2_t arad[],              /* Max alignment adjustment for each object. */
    r2_t step[],              /* Desired adjustment precision for each object. */
    r2_t p[],                 /* (IN/OUT) Corresponding points in each object. */
    double *f2p               /* (OUT) Mismatch for the computed alignment vector. */
  );
  /* Adjusts an alignment vector {p[0..ni-1]} for certain {ni} objects so
    that they are as similar as possible in the neighborhoords of those
    points, as defined by the mismatch function {f2(ni,iscale,p)}.
    Uses iterated quadratic minimization. 
    
    On input, {p[0..ni-1]} must be a guess for the optimum alignment
    vector. On output, {p[0..ni-1]} will be an improved alignment, which
    differs from the given one by at most {arad[i].c[j]}, and differs
    from the optimum by about {step.c[j]}, in each axis {j=0..1}.
    
    If {arad[i].c[j]} is positive, {step[i].c[j]} must be positive. If
    {arad[i].c[j]} is zero or less than {step[i].c[j]}, the coordinate
    {p[i].c[j]} is fixed.
    
    The value of {f2(ni,iscale,p)} is returned on {*f2p}.
    
    The mismatch function {f2} had better have a single minimum within
    the search region, and preferably be approximately quadratic on {p}
    within that region. */

void float_image_align_multi_scale
  ( int ni,                  /* Number of objects to align. */
    float_image_align_mismatch_t *f2, /* Function that evaluates the mismatch between the objects. */
    bool_t quadopt,          /* Use quadratic optimization? */
    r2_t arad[],             /* Max alignment adjustment for each object. */
    r2_t step[],             /* Adjustment step or desired precision for each object. */
    r2_t p[],                /* (IN/OUT) Corresponding points in each object. */
    double *f2p              /* (OUT) Mismatch for the computed alignment vector. */
  );
  /* Adjusts an alignment vector {p[0..ni-1]} for certain {ni} objects
    so that they are as similar as possible in the neighborhoords of
    those points, as defined by the mismatch function {f2(ni,(0,0),p)}. Uses
    a multiscale alignment strategy, with variable step sizes, to
    efficiently search large alignment ranges.
    
    On input, {p[0..ni-1]} must be a guess for the optimum alignment
    vector at scale 0. On output, {p[0..ni-1]} will be an improved
    alignment at scale {(0,0)}, which differs from the given one by at
    most {arad[i].c[j]}, and differs from the optimum by about
    {step.c[j]}, in each axis {j=0..1}.
    
    If {arad[i].c[j]} is positive, {step[i].c[j]} must be positive. If
    {arad[i].c[j]} is zero or less than {step[i].c[j]}, the coordinate
    {p[i].c[j]} is fixed.
    
    The value of {f2(ni,(0,0),p)} is returned on {*f2p}.
    
    The procedure performs an enumerative or quadratic optimization of
    the adjustment vector at various scales {iscale}, starting with
    some coarsest scale {ismax} and ending with scale {(0,0)}.
    Successive scales have at least one of the coordinates reduced by 1.
    
    At each scale, the procedure uses {float_image_align_single_scale_enum}
    or {float_image_align_single_scale_quadopt}, 
    with modified search steps {stp} and modified search radii {rad}
    such that (1) {stp[i].c[j]} is {step[i].c[j]*2^iscale.c[j]}, or
    zero; and (2) {rad[i].c[j]} is about {2*stp[i].c[j]}, or zero. Thus,
    at each scale the adjustment range for {p[i].c[j]} is only a small
    multiple of the adjustment step/precision
    
    The procedure starts the search at a scale {ismax} such that
    {ismax.c[j]} is the smallest integer such that
    {step[i].c[j]*2^iscale.c[j] > arad[i].c[j]}, for all {i} such that
    {arad[i].c[j] > 0}. */

double float_image_align_rel_disp_sqr(int ni, r2_t p[], r2_t q[], r2_t arad[]);
  /* Computes the total squared displacement between {p[0..ni-1]} and 
    {q[0..ni-1]} relative to the radius {arad[i]}, that is,
    
      { SUM { (p[i].c[j] - q[i].c[j])/arad[i].c[j]|^2 : i=0..ni-1,j=0..1 } }
    
    except that terms where {arad[i].c[j]} is zero are ignored. Clients
    may want to add an appropriate multple of this function to the goal
    function in order to bias the search towards the neighborhood of the
    initial guess. */

#endif
