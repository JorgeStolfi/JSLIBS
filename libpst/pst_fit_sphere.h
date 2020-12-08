#ifndef pst_fit_sphere_H
#define pst_fit_sphere_H

/* pst_fit_ellipse.h -- locating a light gauge by fitting an ellipse to its edge. */
/* Last edited on 2009-03-03 14:07:55 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
#include <r2.h>
#include <float_image.h>
#include <ellipse_crs.h>

#include <pst_basic.h>
#include <pst_camera.h>
#include <pst_geom.h>

/* !!! Should adjust the optical center too when given a stretch vector. !!! */

double pst_fit_sphere
  ( float_image_t *IGR, /* Gradient image of a spherical object. */  
    r2_t *Q,      /* Optical center of {IGR}. */
    double *G,    /* (IN/OUT) angular spread of camera. */
    r2_t *K,      /* (IN/OUT) Center of sphere in {IGR}. */
    double *R,    /* (IN/OUT) Radius of sphere. */
    double GAdj,  /* Maximum adjustment allowed in camera spread {G}. */
    double KAdj,  /* Maximum adjustment allowed in {K} coordinates. */
    double RAdj,  /* Maximum adjustment allowed in {R}. */
    int maxIts    /* Max iterations of the optimizer. */
  );
  /* Fits a sphere model to an elliptical outline in a gradient image
    {IGR}.
    
    The procedure assumes that {IGR} is some gradient image derived
    from a photo of a spherical object, which has high values along
    the boundary of the object and low values in the object's interior
    and surround.  (For ordinary photos, the gradient computed by
    {float_image_gradient_relative_sqr} may be adequate for this purpose.)
    
    The procedure also assumes that the center {K} of the sphere is
    anchored and proper; that is, has {K.Z==0}, and its nominal radius
    {R} is less than the focal length {F} of the camera. Since we do
    not model radial distortion yet, the sphere's projection is an
    ellipse. The optical center {Q} of the camera {C} is assumed to
    be fixed as given; and the longest semidiameter is always
    collinear with the vector {K-Q}.
    
    Therefore there are four unknown parameters of the sphere to be
    determined: the Cartesian coordinates {K.X} and {K.Y} of {K}; the
    sphere's nominal radius {R}; and the focal length {F}.
    
    Equivalently, there are four unknown parameters of the elliptical
    projection to be determined: the {X} and {Y} coordinates of its
    center {ctr}, the transverse radius (minor semidiameter) {rad},
    and the length (but not the direction) of the stretch vector
    {E.str}. 
    
    Upon entry, the ellipse {*E} should be an initial guess for the
    sphere's projection, and {*C} should be an initial guess for the
    camera's viewpoint. The initial value of the stretch vector {E.str} is
    ignored, since the stretch is always computed from current values of
    {ctr}, {rad}, and {F}. 
    
    If {KAdj} is nonzero, both coordinates of {K} are adjusted else
    they are both fixed at the given intial values.
    
    If {RAdj} is nonzero, the radius {R} is adjusted, else it is kept
    fixed at the given initial value.
    
    If {GAdj} is nonzero, the camera spead {spr == 1/F} is adjusted,
    else it is maintained at its initial value.
    
    The adjustments to {ctr.X}, {ctr.Y}, and {spr} are absolute, in
    the range {±KAdj}. The adjustments to {rad} are relative (in
    `log scale'), between {rad*t} and {rad/t} where {t =
    (rad+RAdj)/rad}. The same method is used for the adjustments to
    {dst}.
    
    All dimensions and coordinates are in pixels, including the
    focal length {F}. The origin is assumed to be at the
    LOWER left corner of the image, with the Y axis pointing UP.
    
    Uses at most {maxIts} iterations of a non-linear minimizer. The
    goal function is the sum of squared differences between the
    relative gradient of the image and the relative gradient of an
    ideal white ellipse on a black background. */

double pst_fit_sphere_multiscale
  ( float_image_t *IMG, /* Photo of object. */  
    double noise,       /* Std dev of noise in {IMG}, per channel. */
    r2_t *Q,            /* Optical center of {IGR}. */
    double *G,          /* (IN/OUT) angular spread of camera. */
    r2_t *K,            /* (IN/OUT) Center of sphere in {IGR}. */
    double *R,          /* (IN/OUT) Radius of sphere. */
    double GAdj,        /* Maximum adjustment allowed in camera spread {G}. */
    double KAdj,        /* Maximum adjustment allowed in {K} coordinates. */
    double RAdj,        /* Maximum adjustment allowed in {R}. */
    double RMin,        /* Min acceptable radius for multiscale. */
    int maxIts          /* Max iterations of optimizer at initial scale. */
  );
  /* Similar to {pst_fit_sphere}, but uses a multiscale method,
    starting with a scale where the amount of adjustment is between
    {RMin} and {2*RMin}.  Also takes the original 
    image {IMG}, and computes the gradient image {IGR} internally
    with {float_image_gradient_relative_sqr(IMG, noise)}. */

/* INTERNAL TOOLS */

void pst_fit_sphere_data_shrink
  ( double G,
    r2_t K,
    double R, 
    double *Gr,
    r2_t *Kr,
    double *Rr
  );
  /* Returns a camera for a reduced image, given the camera of the original
     image. */

void pst_fit_sphere_data_expand
  ( double G,
    r2_t K,
    double R, 
    double *Gx,
    r2_t *Kx,
    double *Rx
  );
  /* Returns a camera for the original image, given the camera of the
     reduced image. */

#endif
