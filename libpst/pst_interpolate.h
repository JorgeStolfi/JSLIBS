/* pst_interpolate.h -- procedures for interpolating weighted images. */
#ifndef pst_interpolate_H
#define pst_interpolate_H

/* Created on 2011-06-17 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2024-12-22 12:37:34 by stolfi */
/* See the copyright and authorship notice at the end of this file.  */

#include <float_image.h>

/* LINEAR HALFPOINT INTERPOLATION: */

void pst_interpolate_two_values
  ( double v0, double w0,
    double v1, double w1,
    double *vR, double *wR
  );
  /* Estimates a value {*vR} halfway between two points, and its
    reliability weight {*wR}, given the values {v0,v1} at those two
    points, and their reliability weights {w0,w1}. */

void pst_interpolate_two_samples
  ( float_image_t *I, float_image_t *W,
    int32_t c,
    int32_t x0, int32_t y0,
    int32_t x1, int32_t y1,
    double *vR, double *wR
  );
  /* Estimates the value value {*vR} of channel {c} of image {I}
    halfway between the centers of the pixels with indices {x0,y0}
    and {x1,y1}, and its reliability weight {*wR}.
    
    If {W} is given, it must be a single-channel image with same size
    as {I}, containing the reliability weighs of the pixels of {I}.
    If {I} is null, assumes it is all zeros.  If {W} is null,
    assumes it is all ones. */
 
/* CUBIC HALFPOINT INTERPOLATION: */
 
void pst_interpolate_four_values
  ( double vm, double wm,
    double v0, double w0,
    double v1, double w1,
    double vp, double wp,
    double *vR, double *wR
  );
  /* Estimates the value {*vR} at the center of four collinear
    consecutive equidistant points, and its reliability weight {*wR},
    given the values {vm,v0,v1,vp} at those four points (in order) and
    the corresponding reliability weights {wm,w0,w1,wp}. */
  
void pst_interpolate_four_samples
  ( float_image_t *I, float_image_t *W,
    int32_t c,
    int32_t x0, int32_t y0,
    int32_t x1, int32_t y1,
    double *vR, double *wR
  );
  /* Estimates the value {*vR} of channel {c} of image {I}
    halfway between the centers of the pixels with indices {x0,y0}
    and {x1,y1}, and its reliability weight {*wR}. Uses two
    other samples of {I} with indices {xm,ym} and {xp,yp} that
    are collinear with those two points and equally spaced, on both sides.
    
    If {W} is given, it must be a single-channel image with same size
    as {I}, containing the reliability weighs of the pixels of {I}.
    If {I} is null, assumes it is all zeros.  If {W} is null,
    assumes it is all ones. */

#endif
