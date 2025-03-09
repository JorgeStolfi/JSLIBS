#ifndef pst_fit_sphere_H
#define pst_fit_sphere_H

/* pst_fit_sphere.h -- tools for precisely locating spherical gauges in photos. */
/* Last edited on 2025-03-01 19:28:24 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
#include <r2.h>
#include <float_image.h>

#include <pst_basic.h>
#include <pst_light.h>
#include <pst_normal_map.h>
#include <pst_geom.h>

double pst_fit_sphere_and_lights
  ( pst_map_vec_t *IMGV,      /* Photos of a spherical object. */  
    ellipse_crs_t *E,       /* (IN/OUT) Geometric parameters of sphere in {IMGV}. */
    double ctrAdj,          /* Maximum ± adjustment allowed in {ctr} coordinates. */
    double radAdj,          /* Maximum ± adjustment allowed in {rad}. */
    double strAdj,          /* Maximum ± adjustment allowed in {str} coordinates. */
    double adjStep,         /* Ideal adjustment step for geometry parameters. */
    pst_light_vec_t *lhtv,  /* (IN/OUT) light model. */
    int adjustDir,          /* Index of lamp for direction adjustment, or -1. */
    double dirStep,         /* Max lamp direction change per iteration (radians). */ 
    double weightBias,      /* Bias for dark-weighted fitting, or {+INF} for normal fitting. */
    bool_t nonNegative,     /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ,      /* Minimum Z of normal to consider in light fitting. */
    int iterLight,          /* Max geometry-fitting iterations. */
    double tolLight,        /* Convergence criterion for geometry fitting. */
    float_image_t *NRM,     /* (OUT) Normal map of sphere. */ 
    pst_map_vec_t *SYNV       /* (OUT) Synthetic photos of the fitted sphere. */
  );
  /* The procedure assumes that {IMGV[0..NF-1]} (where {NF=IMGV.nel})
    are {NF} monochromatic photos of a spherical object with a
    Lambertian surface of uniform color --- all taken with the same
    camera position and settings, but under different light fields.
    
    The procedure finds the the geometric parameters (center, radius,
    and perspective stretching) and light field parameters (lamp
    direction and intensity) that provide the best fit to those
    photos.
    
    The procedure also assumes that the sphere radius is small
    compared to the camera-to-object distance, so that the perspective
    projection is well approximated by a cylindrical oblique projection.
    Then the sphere's projection is then an ellipse, and the center of the
    sphere projects onto the center {ctr} of the ellipse.
    This assumption should be harmless for a small sphere taken at
    typical camera geometries.

    All dimensions and coordinates are in pixels. The origin is
    assumed to be at the lower left corner of the image, with the
    Y axis pointing up.
    
    INPUTS
    
    Besides the photos {IMGV[0..NF-1]}, the client must provide in the
    variable {E} some initial guesses for the geometric parameters
    of the sphere's projection in those photos.
    
    The client must also set the parameter {lhtv[i].lmpv[j].crad} of
    all lamps in all light fields. That parameter is used in the light
    field fitting loop and in computing the synthetic images. It is
    not changed by this procedure.
    
    OUTPUTS
    
    The procedure will store in {E} the best-fitting sphere
    parameters found. It will also store in {lhtv[i]} the parameters
    of the light field model that provides the best fit between a
    synthetic image of the virtual gauge and the actual photo
    {IMGV[i]}.
    
    The lamp intensities {lhtv[i].src[j].pwr} is computed assuming
    that the object is white (intrinsic lightness is 1.0). If that is
    not the case, the computed value of {srcv[i].pwr} is actually the
    light's intensity multiplied by the object's intrinsic lightness.
    
    Before exiting, the procedure stores in {NRM} the normal map of
    the sphere described by {E}. It also stores in
    {SYNV[i]} the synthetic images of the virtual gauge under 
    the fitted light model {lhtv[i]}.
    
    ALLOCATION

    The client must allocate the vector {lhtv}, including all its
    components, and the output images {NRM}, and {SYNV[0..NF-1]},
    with the proper channel counts and dimensions.
   
    PARAMETER ADJUSTMENT METHOD
    
    The variables {ctrAdj,radAdj,strAdj} specify the maximum amounts
    of adjustment allowed in each parameter. If any of these amounts
    are zero, the corresponding gauge geometry parameters are used
    without any adjustment. If {ctrAdj} is nonzero, both coordinates
    of {E.ctr} are adjusted independently. Ditto for {strAdj} and
    {E.str}.
    
    !!! Revise this, eliminating the {adjStep} parameter: !!!
    
    The {E} parameters are adjusted in equal steps of at most
    {adjStep}. The exact step used may be smaller, in order to ensure
    that all the middle values ({E.rad}, etc.) and extreme values
    ({E.rad ± radAdj}, etc.) of the parameters are sampled.
    
    For each setting of the geometry parameters, the procedure uses an
    iterative least-squares method to adjust the light field
    parameters for each photo, which attempts to minimize the total
    squared difference between the intensities in actual and synthetic
    images. The parameters {adjustDir}, {dirStep}, {weightBias},
    {nonNegative}, {minNormalZ}, {iterLight}, and {tolLight} are used
    in this process; see {pst_fit_light_single_iterative} for details.
    In particular, the light-fitting procedure stops after {iterLight}
    iterations, or when the geometry parameters change by {tolLight}
    or less between two consecutive iterations.
    
    When comparing the actual and synthetic images for light field
    {i}, the procedure ignores any pixel {p} where the intensity of
    the photo {IMGV[i]} is exactly zero. This feature can be used by
    the client to exclude any ``bad'' spots of the input photos from
    the fitting process.
    
    RETURN VALUE
    
    The procedure's return value is the average pixel-by-pixel squared
    intensity difference between the actual photos {IMGV[i]} and the
    synthetic ones {OUTV[i]}, computed over the sphere's
    projection. */

double pst_fit_sphere_evaluate_geometry
  ( pst_map_vec_t *IMGV,      /* Actual photos of a spherical object. */
    ellipse_crs_t *E,       /* Geometric parameters of sphere's projection. */
    pst_light_vec_t *lhtv,  /* (IN/OUT) light model. */
    int adjustDir,          /* Index of lamp for direction adjustment, or -1. */
    double dirStep,         /* Max lamp direction change per iteration (radians). */ 
    double weightBias,      /* Bias for dark-weighted fitting, or {+INF} for normal fitting. */
    bool_t nonNegative,     /* TRUE restricts lamp power and ambient dimming to be non-negative. */
    double minNormalZ,      /* Minimum Z of normal to consider in light fitting. */
    int iterLight,          /* Max iterations for light field fitting. */
    double tolLight,        /* Stopping tolerance for light field fitting. */
    float_image_t *NRM,     /* (OUT) Normal map of best-fit sphere. */
    pst_map_vec_t *SYNV       /* (OUT) Synthetic photos of the fitted sphere. */
  );
  /* Evaluates the parameter set {E}, that describes the presumed
    geometry of a sphere's projection, by comparison with a list of
    actual photos of a spehrical object, {IMGV[0..NF-1]} (where
    {NF=IMGV.nel}). Namely,
  
    * The procedure computes and stores in {NRM} the normal map of the
      sphere with those geometrical parameters.

    * It stores in {SYNV[i]} a sythetic photo obtained from shading
      the normal map {NRM} with the light field {lhtv[i]}, for {i = 0..NF-1}; where
      the the parameters of {lhtv[i].dir} are adjusted so as to
      provide the best match between {SYNV[i]} and {IMGV[i]}.
    
    * It returns the average of sample-by-sample squared intensity
      differences between {IMGV[0..NF-1]} and {SYNV[0..NF-1]} over the
      sphere's projection.
    
    The procedure assumes that {lhtv[i].src[j].crad} is properly set
    for every lamp of every field. The parameters
    {lhtv[i].src[j].dir}, and {lhtv[i].srcv[j].pwr} are estimated by
    least-squares fitting. The domain for least-squares fitting is
    determined iteratively, with at most {iterLight} iterations and
    stopping when the direction changes by {tolLight} or less between
    iterations. The parameters {adjustDir}, {dirStep}, {weightBias},
    {nonNegative}, {minNormalZ} are used in this process; see
    {pst_fit_light_single_iterative} for details. */

double pst_fit_sphere_diff_sqr
  ( float_image_t *NRM,
    float_image_t *AIMG, 
    float_image_t *BIMG
  );
  /* Returns the average of sample-by-sample squared intensity
    differences between {AIMG} and {BIMG}, over the pixels where {NRM}
    is not the null vector. */

#endif
