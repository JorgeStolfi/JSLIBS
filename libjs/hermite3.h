#ifndef hermite3_H
#define hermite3_H

/* Cubic Hermite interpolation */
/* Last edited on 2017-06-25 02:44:05 by stolfilocal */

#define hermite3_H_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

double hermite3_interp(double v0, double d0, double v1, double d1, double a, double b, double t);
  /* Hermite interpolation at some point {t} in the interval {[a _ b]},
    given the values {v0,v1} and the derivatives {d0,d1}
    at the two ends of that interval. */

void hermite3_subsample(int nx, double x[], double dx[], int ns, int ny, double y[]);
  /* Upsamples the sample sequence {x[0..nx-1]} with {ns} steps in each unit interval,
    and stores the result in {y[0..ny-1]}.  Requires {nx>0}, {ns>0}, and 
    {ny == ns*(nx-1)+1}.
    
    Namely, sets {y[i*ns]} to {x[i]} for {i} in {0..nx-1}, and fills the
    other elements {y[k]}, for every {k} in {0..ny-1} that is not
    multiple of {ns}, using a C1 cubic Hermite interpolator. 
    
    If {dx} is not {NULL}, the interpolator assumes that {dx[0..nx-1]} are 
    the derivatives of the function corresponding to the values {x[0..nx-1]}.
    If {dx} is {NULL}, the procedure uses internal estimates of the derivatives
    by {hermite3_estimate_deriv}, that are accurate if the sampled 
    function is any cubic polynomial.  */

void hermite3_estimate_derivs(int nx, double x[], double dx[]);
  /* Estimates derivatives {dx[0..nx-1]} for a function whose
    samples, at equal unit steps, are {x[0..nx-1]}. The 
    estimatives are accurate if the sampled function is any cubic polynomial.
    
    More precisely, if {nx >= 4} the procedure computes the derivative at
    each input sample with {hermite3_estimate_deriv_2_0_2}, except at
    samples {0} and {n-1}, where it uses {hermite3_estimate_deriv_0_1_3}
    and at samples {1} and {n-2}, where it uses {hermite3_estimate_deriv_1_1_2}.
    
    As special cases, it estimates the derivatives by fitting a parabola if {nx == 3},
    an affinity if {nx == 2}, and a constant if {nx == 1}..
    
    Uses {hermite3_estimate_deriv_0_3} for samples {0} and {ny-1},
    {hermite3_estimate_deriv_1_2} for samples {1} and {ny-2},
    and {hermite3_estimate_deriv_2_2} for all the other samples. */ 

double hermite3_estimate_deriv_1_0_1(double vm, double vp);
  /* Assumes that {vm,vp} are two samples of some smooth function,
    taken at arguments {r-1,r+1}, respectively. Estimates the
    derivative at {r} by fitting a quadratic through those points.

    Note that the slope is simply {(vp - vm)/2} and does not
    depend on the function's value at {r}. */ 

double hermite3_estimate_deriv_0_1_2(double v, double vp, double vpp);
  /* Assumes that {v,vp,vpp} are three samples of some smooth function,
    taken at arguments {r,r+1,r+2}, respectively. Estimates the
    derivative at {r} by fitting a quadratic through those points.
    
    Can be used also to estimate the derivative at {r} given values
    {v,vm,vmm} at {r,r-1,r-2}. The result must be negated in that
    case. */ 

double hermite3_estimate_deriv_2_0_2(double vmm, double vm, double vp, double vpp);
  /* Assumes that {vmm,vm,vp,vpp} are four samples of some smooth function,
    taken at arguments {r-2,r-1,r+1,r+2}, respectively. Estimates the
    derivative at {r} by fitting a cubic through those four points. */

double hermite3_estimate_deriv_1_1_2(double vm, double v, double vp, double vpp);
  /* Assumes that {vm,v,vp,vpp} are four samples of some smooth function,
    taken at arguments {r-1,r,r+1,r+2}, respectively. Estimates the
    derivative at {r} by fitting a cubic through the four points.
    
    Can be used also to estimate the derivative at {r} given values
    {vm,v,vp,vpp} at {r+1,r,r-1,r-2}. The result must be negated in that
    case. */

double hermite3_estimate_deriv_0_1_3(double v, double vp, double vpp, double vppp);
  /* Assumes that {v,vp,vpp,vppp} are four samples of some smooth function,
    taken at arguments {r,r+1,r+2,r+3}, respectively. Estimates the
    derivative at {r} by fitting a cubic through the four points.
    
    Can be used also to estimate the derivative at {r} given values
    {v,vp,vpp,vppp} at {r,r-1,r-2,r-3}. The result must be negated in that
    case. */

#endif
