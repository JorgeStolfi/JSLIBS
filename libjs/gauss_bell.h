#ifndef gauss_bell_H
#define gauss_bell_H

/* Tools for parameters of a Gaussian distribution. */
/* Last edited on 2019-12-05 10:24:27 by jstolfi */

#include <bool.h>

typedef struct gauss_bell_parms_t 
  { double avg; /* Mean value.*/
    double dev; /* Standard deviation. */
  } gauss_bell_parms_t;
  /* Parameters of a generic Gaussian bell function or distribution,
     with mean {avg} and standard deviation {dev}. */

#define gauss_bell_HUGE_ARG (38.6)
  /* A value such that {|z| >= gauss_bell_HUGE_ARG} implies
    {exp(-z^2/2)} underflows. */

#define gauss_bell_BIG_ARG (8.6)
  /* A value such that {|z| >= gauss_bell_BIG_ARG} implies
    {exp(-z^2/2) <= 1.0e-16} . */
 
#define gauss_bell_TINY_ARG (1.4e-8)
  /* A value such that {|z| <= gauss_bell_TINY_ARG} implies
    {exp(-z^2/2) >= 1 - 1.0e-16}. */
 
double gauss_bell_eval(double x, double avg, double dev);
  /* The UNNORMALIZED Gaussian bell {exp(-((x - avg)/dev)^2/2)}, 
     but safer/faster for big and small args.  
     
     The sign of {dev} is ignored.  If {dev} is zero, 
     returns 1 if {x=avg}, 0 otherwise.  If {dev} is {±INF},
     the result is 1.0 for any {z}. */

#endif
