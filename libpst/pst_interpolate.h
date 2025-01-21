/* pst_interpolate.h -- procedures for interpolating weighted images. */
#ifndef pst_interpolate_H
#define pst_interpolate_H

/* Created on 2011-06-17 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-01-16 12:54:52 by stolfi */
/* See the copyright and authorship notice at the end of this file.  */

#include <float_image.h>

/* LINEAR HALFPOINT INTERPOLATION: */

void pst_interpolate_two_values
  ( double v0, double w0,
    double v1, double w1,
    double *vRP, double *wRP
  );
  /* Estimates a value {vR} halfway between two points, and its
    reliability weight {wR}, given the values {v0,v1} at those two
    points, and their reliability weights {w0,w1}. Returns the 
    results in {*vRP} and {*wRP}. */

/* CUBIC HALFPOINT INTERPOLATION: */
 
void pst_interpolate_four_values
  ( double vm, double wm,
    double v0, double w0,
    double v1, double w1,
    double vp, double wp,
    double *vRP, double *wRP
  );
  /* Estimates the value {vR} at the center of four collinear
    consecutive equidistant points, and its reliability weight {wR},
    given the values {vm,v0,v1,vp} at those four points (in order) and
    the corresponding reliability weights {wm,w0,w1,wp}. Returns the 
    results in {*vRP} and {*wRP}. */

#endif
