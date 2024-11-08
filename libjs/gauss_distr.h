#ifndef gauss_distr_H
#define gauss_distr_H

/* Tools for parameters of a Gaussian distribution (with integral 1.0). */
/* Last edited on 2024-11-05 20:56:37 by stolfi */

#include <bool.h>
#include <gauss_bell.h>

typedef gauss_bell_parms_t gauss_distr_parms_t;
  /* Parameters of a generic Gaussian distribution, namely 
     mean {.avg} and standard deviation {.dev}. */
     
#define gauss_distr_HUGE_ARG (gauss_bell_HUGE_ARG)
  /* A value such that {|z| >= gauss_distr_HUGE_ARG} implies 
    that the standard PDF at {z} underflows. */
     
#define gauss_distr_BIG_ARG (gauss_bell_BIG_ARG)
  /* A value such that {|z| >= gauss_distr_BIG_ARG} implies 
    that the standard PDF at {z} is less than {10^{-16}}. */
 
#define gauss_distr_TINY_ARG (gauss_bell_TINY_ARG)
  /* A value such that {|z| <= gauss_distr_TINY_ARG} implies
    that the standard PDF at {z} is greater than {1 - 10^{-16}}. */

double gauss_distr_PDF(double x, double avg, double dev);
  /* The probability density function of a normal deviate with 
    mean {avg} and deviation {dev}, evaluated at {x}.
    
    Namely, returns {A*exp[-((x-avg)/dev)^2/2]} where {A = 1/(|dev|*sqrt(2*PI))}
    is such that the integral is over all {x} is 1.0.  Possibly safer and/or
    faster for very large or very small {x}.  
    
    The sign of {dev} is ignored. If {dev} is {±INF},  
    returns 0 for any {x}. If {dev} is zero, returns {+INF} 
    for {x=avg} (or {NAN} if {avg} is {±INF}) and 0 otherwise.  */
    
double gauss_distr_CDF(double x, double avg, double dev);
  /* The cumulative distribution function of a Gaussian deviate 
    with average {avg} and deviation {dev}.  Namely, the integral
    of {gauss_distr_PDF(y, avg, dev)} for {y} from {-INF} to {x}.
    
    The sign of {dev} is ignored. If {dev} is {±INF},
    returns 0 {x < +INF} and 1 for {x=+INF}.  If {dev} is zero,
    returns 0 for  {x<avg}, 1 for {x>avg}, and {NAN} for {x=avg}. */
    
double gauss_distr_integral(double x0, double x1, double avg, double dev);
  /* Returns the integral from {x0} to {x1} of the normalized Gaussian
    PDF with average {avg} and deviation {dev}. Namely, the probability
    that a Gaussian deviate with those parameters will be between {x0}
    and {x1}.
     
    The sign of {dev} is ignored. 
    
    If {dev} is zero, returns {+1} if {x0<avg<x1}, {-1} if {x1<avg<x0},
    0 if {x0,x1} are on the same side of avg, and {NAN} if either of
    them is equal to {avg}.
    
    If {dev} is infinte, returns {+1} if {x0=-INF} and {x1=+INF}, {-1}
    if {x0=+INF} and {x1=-INF}, {NAN} if only one of them is infinite,
    and zero otherwise. 
    
    If {dev} is finite and nonzero, and {x0=x1}, returns 0. */

#endif
