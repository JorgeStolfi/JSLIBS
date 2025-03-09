/* pst_interpolate.h -- procedures for interpolating weighted images. */
#ifndef pst_interpolate_H
#define pst_interpolate_H

/* Created on 2011-06-17 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-02-20 08:26:20 by stolfi */
/* See the copyright and authorship notice at the end of this file.  */

#include <float_image.h>

/* LINEAR HALFPOINT INTERPOLATION: */
 
void pst_interpolate_four_values
  ( double vm, double wm,
    double v0, double w0,
    double v1, double w1,
    double vp, double wp,
    double *vR_P, double *wR_P
  );
  /* Estimates the value {vR} at the center of four collinear
    consecutive equidistant points, and its reliability weight {wR},
    given the values {vm,v0,v1,vp} at those four points (in order) and
    the corresponding reliability weights {wm,w0,w1,wp}. Returns the 
    results in {*vR_P} and {*wR_P}.
    
    The weights must be finite non-negative numbers. If any weight is
    nonzero, the corresponding value must be finite. If {w0} and {w1}
    are both positive, uses linear, quadratic, or cubic interpolation,
    depending on whether none, one, or two of {wm,wp} are positive.
    Otherwise the result {vR} will be {NAN} and its weight {wR} will be
    zero. Idem if overflow or underflow occurs. */

#endif
