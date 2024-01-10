/* Routines for standard interval arithmetic    */
/* Created by Jorge Stolfi 93-01-13             */

#ifndef IA_H
#define IA_H

#include <flt.h>
#include <stdio.h>
#include <math.h>

typedef struct {Float lo, hi; } Interval;

#define IA_FULL  (ia_full())
#define IA_ISFULL(x)  ((x.lo <= MinusInfinity) || (x.hi >= PlusInfinity))
#define IA_NORMFULL(x)  if (IA_ISFULL(x)) x = IA_FULL

Interval ia_full (void); /* [-Inf .. +Inf] */

Interval ia_const(Float x, Float err); /* $x$ plus or minus $err$ */

Interval ia_add   (Interval x, Interval y);
Interval ia_sub   (Interval x, Interval y);
Interval ia_neg   (Interval x);
Interval ia_scale (Interval x, Float alpha, Float zeta); /* $\alpha x / \zeta$ */
Interval ia_shift (Interval x, Float gamma);  /* $x + \gamma$ */
Interval ia_mul   (Interval x, Interval y);
Interval ia_div   (Interval x, Interval y);
Interval ia_sqr   (Interval x); /* same as ia_mul(x,x), only better */
Interval ia_inv   (Interval x);
Interval ia_sqrt  (Interval x);

Interval ia_affine (
    Interval x, Float alpha, 
    Float zeta,
    Float gamma, Float delta
  );
  /* 
    Computes $\alpha x / \zeta + \gamma \pm \delta. */

Interval ia_affine_2(
    Interval x, Float alpha,
    Interval y, Float beta, 
    Float zeta, 
    Float gamma, Float delta
  );
  /* 
    Computes $(\alpha x + \beta y)/ \zeta + \gamma \pm \delta. */

Interval ia_meet (Interval x, Interval y);
  /* 
    Intersection of $x$ and $y$; error if disjoint. */

void ia_print (FILE *f, Interval x);

#endif
