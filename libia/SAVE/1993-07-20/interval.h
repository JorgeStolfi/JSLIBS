/* Routines for ordinary interval arithmetic */

#ifndef INTERVAL_H
#define INTERVAL_H

#include "foifloat.h"
#include <stdio.h>
#include <math.h>

typedef struct {Float lo, hi; } Interval;

#define IV_FULL  (iv_full())
#define IV_ISFULL(x)  ((x.lo <= MinusInfinity) || (x.hi >= PlusInfinity))
#define IV_NORMFULL(x)  if (IV_ISFULL(x)) x = IV_FULL

Interval iv_full (void); /* [-Inf .. +Inf] */

Interval iv_const(Float x, Float err); /* $x$ plus or minus $err$ */

Interval iv_add (Interval x, Interval y);
Interval iv_neg (Interval x);
Interval iv_sub (Interval x, Interval y);
Interval iv_mul (Interval x, Interval y);
Interval iv_sqr (Interval x); /* same as iv_mul(x,x), only better */
Interval iv_inv (Interval x);
Interval iv_sqrt (Interval x);

Interval iv_affine (
    Interval x,
    Float alpha, Float beta, Float gamma
  );
  /* Computes $x*\alpha + \beta \pm \gamma$. */

Interval iv_meet (Interval x, Interval y);
  /* Intersection of $x$ and $y$; error if disjoint. */

void iv_print (FILE *f, Interval x);

#endif
