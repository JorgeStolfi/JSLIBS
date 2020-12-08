#ifndef pst_fit_ellipse_H
#define pst_fit_ellipse_H

/* pst_fit_ellipse.h -- locating a light gauge by fitting an ellipse to its edge. */
/* Last edited on 2009-08-22 13:04:17 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
#include <r2.h>
#include <ellipse_crs.h>
#include <float_image.h>

#include <pst_basic.h>
#include <pst_light.h>
#include <pst_normal_map.h>
#include <pst_geom.h>

/* !!! Split off the gradient computation !!! */

double pst_fit_ellipse
  ( float_image_t *IGR, /* Gradient image of a spherical object. */  
    ellipse_crs_t *EP,  /* (IN/OUT) Geometric parameters of sphere in {IGR}. */
    double ctrAdj,      /* Maximum adjustment allowed in {EP.ctr} coordinates. */
    double radAdj,      /* Maximum adjustment allowed in {EP.rad}. */
    double strAdj,      /* Maximum adjustment allowed in {EP.str} coordinates. */
    int maxIts          /* Max iterations of the optimizer. */
  );
  /* Fits an ellipse to the outline of a spherical object in a
    gradient image {IGR}.
    
    The procedure assumes that {IGR} is some gradient image derived
    from a photo of a spherical object, which has high values along
    the boundary of the object and low values in the object's interior
    and surround.  (For ordinary photos, the gradient computed by
    {float_image_gradient_relative_sqr} may be adequate for this purpose.)
    
    The procedure also assumes that the object's projection is an
    ellipse.  This assumption is correct even for conical projection
    and is practically correct even if the camera has some radial
    distortion. 
    
    Upon entry, the procedure assumes that {*EP} is an initial guess
    for the geometric parameters (center, radius, and perspective
    stretching) of the shpere's outline in the photo. The adjusted
    parameters of the ellipse are returned in {*EP}.

    If {ctrAdj} is nonzero, both coordinates of the center are adjusted,
    else they are both fixed at the given intial values.
    
    If {radAdj} is nonzero, the radius is adjusted, else
    it is kept fixed at the given initial value.
    
    If {strAdj} is nonzero, both coordinates of the stretch vector,
    are adjusted, else both are fixed at the given initial values.
    
    The adjustments to {EP.ctr} and {EP.str} are absolute, in the
    range {±ctrAdj} and {±strAdj}, respectively. The adjustments
    to {rad} are relative (in `log scale'), between {EP.rad*t} and
    {EP.rad/t} where {t = (EP.rad+radAdj)/rad}.

    All dimensions and coordinates are in pixels. The origin is
    assumed to be at the corner of the image's domain 
    nearest to pixel [0,0], with both axes pointing towards
    increasing indices.
    
    Uses at most {maxIts} iterations of a non-linear minimizer. The
    goal function is the sum of squared differences between the
    relative gradient of the image and the relative gradient of an
    ideal white ellipse on a black background. */

double pst_fit_ellipse_multiscale
  ( float_image_t *IMG, /* Photo of object. */  
    double noise,       /* Std dev of noise in {IMG}, per channel. */
    ellipse_crs_t *EP,  /* (IN/OUT) Geometric parameters of sphere in {IMG}. */
    double ctrAdj,      /* Maximum adjustment allowed in {EP.ctr} coordinates. */
    double radAdj,      /* Maximum adjustment allowed in {EP.rad}. */
    double strAdj,      /* Maximum adjustment allowed in {EP.str} coordinates. */
    double minRadius,   /* Min acceptable radius for multiscale. */
    int maxIts          /* Max iterations of optimizer at initial scale. */
  );
  /* Similar to {pst_fit_ellipse}, but uses a multiscale method,
    starting with a scale where the amount of adjustment is between
    {minRadius} and {2*minRadius}.  Also takes the original 
    image {IMG}, and computes the gradient image {IGR} internally
    with {float_image_gradient_relative_sqr(IMG, noise)}.
    
    The iteration limit {maxIts} applies to the optimization at the 
    coarser scale only.  For the intermediate scales, the limit 
    is either {maxIts} or a fixed internal limit, whichever is 
    smaller.  Thus, if {maxIts} is zero, the optimization 
    will be supressed at all scales. */

double pst_fit_ellipse_eval(float_image_t *IGR, ellipse_crs_t *EP);
  /* A quadratic-like function that expresses the mismatch between
    the normalized gradient image {IMG} and the ideal ellipse
    whose geometry is {*EP}. */

/* Internal tools: */

#define pst_fit_ellipse_nw 3
  /* Size of filter window for {pst_fit_ellipse_image_shrink}. */

float_image_t *pst_fit_ellipse_image_shrink(float_image_t *IMG);
  /* Reduces the image {IMG} for multiscale fit. */

ellipse_crs_t pst_fit_ellipse_geom_shrink(ellipse_crs_t *EP);
  /* Produces the specifications {grd} of an ellipse that results from
    the ellipse described by {EP} by the same domain reduction of
    {pst_fit_ellipse_image_shrink}. */

ellipse_crs_t pst_fit_ellipse_geom_expand(ellipse_crs_t *EP);
  /* Produces the specifications {gex} of an ellipse that results from
    the ellipse described by {EP} by the inverse of the domain reduction of
    {pst_fit_ellipse_image_shrink}. */

#endif
