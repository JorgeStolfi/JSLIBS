#ifndef interp_spline_H
#define interp_spline_H

/* Spline smoothing/interpolation of equally-spaced numerical series. */
/* Last edited on 2013-10-25 23:25:39 by stolfilocal */ 

#include <bool.h>
#include <ix.h>

/* The procedures in this interface convert a sequence of real samples
  {s[0..ns-1]} into a real-valued function {f} of a real argument {z},
  by some local smoothing or interpolation method.
  
  In general, the value {f(z)} is some weighted average of a certain
  number {nw} of consecutive data samples {s[k..k+nw-1]} with certain
  weights {w[0..nw-1]}; where
  
    * the number {nw} (the /number of taps/ of the method) 
    depends only on two discrete parameters, the /kind/ {knd} and 
    the /order/ {ord};
   
    * the initial index {k} is such that {z} is in the range 
    {k+(nw/2)+[-1/2 _ +1/2]}, that is, {k = floor(z - (nw-1)/2)};
    
    * the weights {w[0..nw-1]} depend only on {knd}, {ord}, and the 
    fractional part of {z}.
  
  These rules apply only when all the indices {k..k+nw-1} lie in the 
  range {0..ns-1}.  In general, each `raw' sample index {k+i} is mapped 
  to a `reduced' index {ix[i] = ix_reduce(k+i,ns,red)}, for a specific
  index reduction method {red} (see {indexing.h}). If {ix[i]} is {-1}
  the corresponding data sample is omitted from the average.
  Therefore, the interpolating function is given by
  
    { f(z) = SUM{w[i]*s[ix[i]]} / SUM{w[i]} }
    
  where both sums are taken over the values {i} in {0..nw-1} such that
  {ix[i]} is non-negative.  Note that the value is {NAN} if
  all indices {ix[0..nw-1]} are {-1}, for example when
  {red = ix_reduction_SINGLE} and {z} is sufficiently far from
  the interval {[0 _ ns]}.
    
  Note that, for any integer {d}, the computation of{f(z+d)} uses the
  same weights {w[0..nw-1]} as {f(z)}, but applied to the samples with
  raw indices {k+d..k+d+nw-1} instead of {k..k+nw-1}. Therefore, when
  downsampling {s} with a fixed integer number {m} of regularly spaced
  sub-samples per data sample, one can precompute the {m} sets of
  distinct weights needed and reuse them {ns} times. */

typedef enum 
  { interp_spline_kind_B,    /* B-splines (convex, broad window, not interpolating). */
    interp_spline_kind_I,    /* I-splines (interpolating, min window, not convex). */
    interp_spline_kind_O,    /* O-splines (interpolating, min degree, not convex). */
    interp_spline_kind_NUM   /* Number of spline kinds. */
  } interp_spline_kind_t;

int interp_spline_compute_num_samples(int ord, interp_spline_kind_t knd);
  /* Computes the number of consecutive data samples {nw} needed for
    interpolation at a generic real argument specified by {ord} and {knd}. */

void interp_spline_get_indices(double z, int ns, ix_reduction_t red, int nw, int ix[]);
  /* Computes the indices {ix[0..nw-1]} of the {nw} samples needed to perform
    interpolation in a sequence of samples at the position {z},
    assuming that the nominal position of a sample with index {k} is {z = k+0.5}.
    
    Indices {ix[k]} that fall outside the range {0..ns-1} are reduced
    according to {ix_reduce(ix, ns, red)}. Some reduced indices {ix[k]} may
    be {-1} to denote `no such sample'. */

void interp_spline_get_weights(double z, int ord, interp_spline_kind_t knd, int nw, double w[]);
  /* Computes the sample weights {w[0..nw-1]} for evaluation of {f(z)}.
    Assumes that {ord>=-1} and {nw=interp_spline_compute_num_samples(ord, kind)}. */

#endif
