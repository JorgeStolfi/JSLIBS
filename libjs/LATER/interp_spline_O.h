#ifndef ospline_interp_H
#define ospline_interp_H

/* Spline interpolation of equally-spaced numerical series. */
/* Last edited on 2024-11-18 11:19:47 by stolfi */ 

#include <stdint.h>

#include <bool.h>

/* This module defines an interpolation method, depending on a paramter
  {ord >= -1}, that converts a sequence of samples {s[0..ns-1]} into a
  spline (piecewise polynomial function) {f} of degree {ord+1} which has
  continous derivatives of any order in {0..ord}.

  Assuming that each data sample with index {i} is located at{z=i+0.5},
  the spline interpolates the given samples; that is, {f(i+0.5) = s[i]}
  for {i} in {0..ns-1}.
  
  In particular, when {ord} is {-1}, {f} is the nearest-sample
  interpolating staircase. When {ord} is 0, {f} is the linear
  interpolation of consecutive samples. */
  
uint32_t interp_spline_O_compute_num_samples(int32_t ord);
  /* Computes the number {nw} of consecutive samples needed for
    interpolation of order {ord} at a general point {z}. */

void interp_spline_O_get_weights(double z, int32_t ord, uint32_t nw, double wt[]);
  /* Computes the indices {i[0..nw-1]} for evaluation of {f(z)}.
    Assumes that {ord >= -1} and {nw} is 
    {interp_spline_O_compute_num_samples(ord)}. */

#endif
