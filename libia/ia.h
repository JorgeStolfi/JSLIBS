/* Routines for standard interval arithmetic    */
/* Last edited on 2024-12-31 00:48:01 by stolfi */
/* Created by Jorge Stolfi 93-01-13             */

#ifndef ia_H
#define ia_H

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include <flt.h>
#include <jsrandom.h>

typedef struct {Float lo, hi; } Interval;

void ia_init (void);
  /* Initializes the constants (ia_full, etc.) */
  /* MUST be called at least once before any of the routines below. */

Interval ia_full (void); /* The "anything" interval, [-Inf .. +Inf] */

int32_t ia_is_full (Interval *x); 
  /* True iff {*x.lo} is -Inf or {*x.hi} is +Inf */
  
void ia_norm_full (Interval *x);
  /* If {*x.lo} is -Inf or {*x.hi} is +Inf, sets {*x} to [-Inf .. +Inf] */

Interval ia_const (Float x, Float err); /* {x} plus or minus {err} */

Interval ia_int_const (int32_t i);  /* {i}, with rounding error if too big. */

Interval ia_add   (Interval x, Interval y);
Interval ia_sub   (Interval x, Interval y);
Interval ia_neg   (Interval x);
Interval ia_scale (Interval x, Float alpha, Float zeta); /* {alpha * x / zeta} */
Interval ia_shift (Interval x, Float gamma);  /* {x + gamma} */
Interval ia_mul   (Interval x, Interval y);
Interval ia_div   (Interval x, Interval y);
Interval ia_sqr   (Interval x); /* same as {ia_mul(x,x)}, only better */
Interval ia_inv   (Interval x);
Interval ia_exp   (Interval x);
Interval ia_sqrt  (Interval x);
Interval ia_abs   (Interval x);
Interval ia_max   (Interval x, Interval y);
Interval ia_min   (Interval x, Interval y);

Interval ia_affine (
    Interval x, Float alpha, 
    Float zeta,
    Float gamma, Float delta
  );
  /* Computes {alpha * x / zeta + gamma ± delta}. */

Interval ia_affine_2(
    Interval x, Float alpha,
    Interval y, Float beta, 
    Float zeta, 
    Float gamma, Float delta
  );
  /* Computes {(alpha * x + beta * y)/ zeta + gamma ± delta}. */

/*** MISCELLANEOUS TOOLS ***/

Float ia_mid (Interval x);
  /* Approximate midpoint of {x}, guaranteed to be in {x}.
    Finite as long as {x} is finite. */

Float ia_rad (Interval x);
  /* Radius of {x} from its midpoint {m = ia_mid(x)}.
    I.e. a value {r} such that {[m-r _ m+r]} contains {x}.
    Finite as long as {x} is finite. */

Interval ia_meet (Interval x, Interval y);
  /* Intersection of {x} and {y}; error if disjoint. */

Interval ia_join (Interval x, Interval y);
  /* Smallest interval contining {x} and {y}. */

Interval ia_interp (Float x0, Interval y0, Float x1, Interval y1, Float x);
  /* Returns an interval {y} that contains {r(x)} for any affine function
    (straight line) such that {r(x0)} is in {yp} and {r(x1)} is in {y1}. */

Interval ia_throw (void);
  /* Returns a random interval, suitable for testing.
    The client must have called {srandom(<seed>)}. */

void ia_print (FILE *wr, Interval x);
  /* Prints {x} to file {wr}, in the format "[ {x.lo} __ {x.hi} ]". */
    
void ia_print_bound(FILE *wr, Float v, int32_t which, int32_t full);
  /* If {full} is 0, prints the value {v} to file {wr}; 
    rounded down if {which} is 0, rounded up if {which} is 1.
    If {full} is 1, prints a bunch of '*'s instead. */

#endif
