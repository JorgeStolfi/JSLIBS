#ifndef ospline_interp_H
#define ospline_interp_H

/* Spline interpolation of equally-spaced numerical series. */
/* Last edited on 2013-10-25 23:22:14 by stolfilocal */ 

#include <bool.h>

/* This module defines an interpolation method, depending on a paramter
  {ord >= -1}, that converts a sequence of samples {s[0..ns-1]} into a
  spline (piecewise polynomial function) {f} of degree {ord+1} which has
  continous derivatives of any order in {0..ord}.

  Assuming that each data sample with index {i} is located at {z =
  i+0.5}, then the spline interpolates the given samples;
  that is, {f(i+0.5) = s[i]} for {i} in {0..ns-1}.
  
  In particular, when {ord} is {-1}, {f} is the nearest-sample interpolating staircase.
  When {ord} is 0, {f} is the linear interpolation of consecutive samples.  */
  
int interp_spline_O_compute_num_samples(int ord);
  /* Computes the number of samples {nw} needed along each axis for
    interpolation of order {ord}. */

void interp_spline_O_get_weights(double z, int ord, int nw, double wt[]);
  /* Computes the indices {i[0..nw-1]} for evaluation of {f(z)}.
    Assumes that {ord >= -1} and {nw} is 
    {interp_spline_O_compute_num_samples(ord)}. */

#endif
