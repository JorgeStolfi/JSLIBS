#ifndef jsroot_H
#define jsroot_H

#include <stdint.h>

/* Routines to find roots of univariate functions. */
/* Last edited on 2008-10-05 18:05:06 by stolfi */

double jsroot_bisect_secant(double xLo, double xHi, double (*f)(double x));
  /* Finds a {x} such that {f(x)} is approximately 0, by alternate
    bisection and secant steps. Assumes that {f} is continuous.
    Requires that {f(xLo) * f(xHi) <= 0}. */


#endif
