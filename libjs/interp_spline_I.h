#ifndef ispline_interp_H
#define ispline_interp_H

/* Minimum-width spline interpolation of equally-spaced numerical series. */
/* Last edited on 2024-11-18 11:22:15 by stolfi */ 

#include <stdint.h>

#include <bool.h>

/* This module defines an interpolation method, depending on a paramter
  {ord >= -1}, that converts a sequence of samples {s[0..ns-1]} into a
  spline (piecewise polynomial function) {f} of degree {ord+1} which has
  continous derivatives of any order in {0..ord}.

  Assuming that each data sample with index {i} is located at {z =
  i+0.5}, then the spline interpolates the given samples;
  that is, {f(i+0.5) = s[i]} for {i} in {0..ns-1}.
  
  When {ord} is {-1}, {f} is the nearest-sample interpolating staircase.
  
  When {ord} is 0, {f} is the linear interpolation of consecutive samples.
  
  When {ord} is 1, {f} is formed by quadratic arcs, each spanning an interval
  {[i _ i+1/2]} or {[i-1/2 _ i]} for some integer {i}. 
  
  These functions are not yet implemented for {ord>=2}.*/
  
uint32_t interp_spline_I_compute_num_samples(int32_t ord);
  /* Computes the number of samples {m} needed for interpolation of the
    specified order {ord}. Assumes {ord >= -1}. */

void interp_spline_I_get_weights(double z, int32_t ord, uint32_t nw, double wt[]);
  /* Computes the weights {wt[0..nw-1]} for evaluation of {f(z)}.
    Assumes that {ord >= -1} and {nw} is
    {interp_spline_I_compute_num_samples(ord)}. */

#endif
