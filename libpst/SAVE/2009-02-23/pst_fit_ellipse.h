#ifndef pst_fit_ellipse_H
#define pst_fit_ellipse_H

/* pst_fit_ellipse.h -- locating a light gauge by fitting an ellipse to its edge. */
/* Last edited on 2009-02-22 14:34:32 by stolfi */

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

double pst_fit_ellipse
  ( float_image_t *IMG, /* Monochromatic photo of object. */  
    ellipse_crs_t *EP,  /* (IN/OUT) Geometric parameters of sphere in {IMGV}. */
    double ctrAdj,      /* Maximum ± adjustment allowed in {ctr} coordinates. */
    double radAdj,      /* Maximum ± adjustment allowed in {rad}. */
    double strAdj,      /* Maximum ± adjustment allowed in {str} coordinates. */
    int maxIts          /* Max iterations of the optimizer. */
  );
  /* The procedure assumes that {IMG} is a monochromatic photo of a
    spherical object against a contrating background. Upon entry, the
    procedure assumes that {*EP} is an initial guess for the geometric
    parameters (center, radius, and perspective stretching) of the
    shpere's outline in the photo. The procedure adjusts {*EP} so that
    it best fits the photo.
    
    The procedure assumes that the sphere radius is small compared to
    the camera-to-object distance, so that the perspective projection
    is well approximated by a parallel oblique projection. Then the
    sphere's projection is an ellipse, and the center of the sphere
    projects onto the center {ctr} of the ellipse. This assumption
    should be harmless for a small sphere taken at typical camera
    geometries.
    
    If {ctrAdj} is nonzero, both coordinates of the center are adjusted,
    else they are both fixed at the given intial values.
    
    If {radAdj} is nonzero, the radius is adjusted, else
    it is kept fixed at the given initial value.
    
    If {strAdj} is nonzero, both coordinates of the stretch vector,
    are adjusted, else both are fixed at the given initial values.

    All dimensions and coordinates are in pixels. The origin is
    assumed to be at the lower left corner of the image, with the
    Y axis pointing up.
    
    Uses at most {maxIts} iterations of a non-linear minimizer. The
    goal function is the sum of squared differences between the
    relative gradient of the image and the relative gradient of an
    ideal white ellipse on a black background. */

double pst_fit_ellipse_multiscale
  ( float_image_t *IMG,     /* Monochromatic photo of object. */  
    ellipse_crs_t *EP,      /* (IN/OUT) Geometric parameters of sphere in {IMGV}. */
    double ctrAdj,          /* Maximum ± adjustment allowed in {ctr} coordinates. */
    double radAdj,          /* Maximum ± adjustment allowed in {rad}. */
    double strAdj,          /* Maximum ± adjustment allowed in {str} coordinates. */
    double minRadius,       /* Min acceptable radius for multiscale. */
    int maxIts              /* Max iterations of optimizer at initial scale. */
  );
  /* Similar to {pst_fit_ellipse}, but uses a multiscale method,
    starting with a scale where the amount of adjustment is between
    {minRadius} and {2*minRadius}. */

double pst_fit_ellipse_eval(float_image_t *IGR, ellipse_crs_t *EP);
  /* A quadratic-like function that expresses the mismatch between
    the normalized gradient image {IMG} and the ideal ellipse
    whose geometry is {*EP}. */

#endif
