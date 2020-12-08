#ifndef jsqroots_H
#define jsqroots_H

#include <stdint.h>
#include <jsmath.h>

/* A quadratic equation solver. */
/* Last edited on 2008-07-16 01:26:45 by stolfi */
    
int32_t roots_quadratic(double A, double B, double C, double *r1, double *r2, double *im);
  /* Finds the roots of a quadratic equation {A*x*x + B*x + C == 0}.
    If the roots are real, stores them in {*r1} and {*r2} (so that
    {*r1 <= *r2}) and sets {*im} to 0. If they are complex, sets {*r1}
    and {*r2} to the real part, and {*im} to the (non-negative)
    imaginary part. Returns the sign ({-1}, {0}, or {+1}) of the
    determinant as the result.
    
    The function tries to return sensible results even in degenerate
    cases. If {A} is zero and {B} is nonzero, sets {*r1} to the root
    of the 1st degree equation {B*x +C == 0}, {*r2=NaN}, and {*im=0}.
    If {A} and {B} are both zero, sets {*r1=*r2=NaN}, {*im=0}. If
    any coefficient is infinte or {NaN}, sets {*r1=*r2=im=NaN}. */

int32_t roots_proper_quadratic(double A, double B, double C, double *r1, double *r2, double *im);
 /* Like {roots_quadratic}, but fails if any of {A,B,C} is zero. */

/* Created by J. Stolfi, IC-UNICAMP on 2008-07-14 */

#endif
