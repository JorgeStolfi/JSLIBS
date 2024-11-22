#ifndef bspline_interp_H
#define bspline_interp_H

/* B-spline smoothing of equally-spaced numerical series. */
/* Last edited on 2024-11-18 09:19:26 by stolfi */ 

#include <stdint.h>

#include <bool.h>

/* B-spline smoothing of a given order {ord >= -1} converts a sequence of
  samples {s[0..ns-1]} into a spline (piecewise polynomial function) {f}
  of degree {ord+1} which has continous derivatives of any order in
  {0..ord}.
  
  Each piece of the spline is a polynomial of degree {ord+1} defined on
  a unit-length interval; specifically, {[i _ i+1]} if {ord} is odd,
  {[i-1/2 _ i+1/2]} if ord is even, for some integer {i}.
  
  Assuming that each data sample with index {i} is located at {z =
  i+0.5}, then the joints of the spline are located halfway between
  consecutive samples if {ord} is odd, at the samples if {ord} is even.

  The spline {f} is interpolatory (that is, {f(i+0.5) = s[i]}) when
  {ord} is {-1} (when it is the nearest-sample interpolating staircase)
  or {0} (when it is the linear interpolation of consecutive samples).
  For {ord>=1}, the spline {f} is interpolatory if the samples come from
  a polynomial of degree at most {ord+1}, and in that case {f} is that
  polynomial. Otherwise the spline is not interpolatory in general. In
  any case, the spline value is a convex combination of the {nw} samples
  used in the interpolation. */

uint32_t interp_spline_B_compute_num_samples(int32_t ord);
  /* Computes the number of samples {nw} needed for B-spline
    interpolation with continuity {ord} (namely, {ord+2}).
    Assumes {ord >= -1}. */

void interp_spline_B_get_weights(double z, int32_t ord, uint32_t nw, double wt[]);
  /* Computes the sample weights {wt[0..nw-1]} for evaluation of {f(z)}.
    Assumes that {ord >= -1} and {nw} is 
    {interp_spline_B_compute_num_samples(ord)}. */

#endif
