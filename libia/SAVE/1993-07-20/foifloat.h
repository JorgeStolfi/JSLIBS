/* basic floating-point definitions */

#ifndef FOIFLOAT_H
#define FOIFLOAT_H

#include <values.h>
#include <floatingpoint.h>
#include <math.h>
#include <stdio.h>

/* Type of Interval bounds, FOI coefficients, etc.: */
typedef
  float Float;

typedef
  float *FloatP;

/* Default $printf$ parameters for type Float: */

#define F_FMT_DIGITS 8
#define F_FMT_WIDTH (F_FMT_DIGITS + 6)
#define F_FMT_DECIMALS (F_FMT_DIGITS - 1)
#define F_FMT_SPEC "%14.7e"

/* Common constants for type Float: */

#define One     (1.0)
#define Zero    (0.0)
#define Half    (0.5)
#define Quarter (0.25)
#define Two     (2.0)
#define Three   (3.0)
#define Four    (4.0)
#define MaxFloat MAXFLOAT
#define Infinity  HUGE_VAL
#define MinusInfinity (-Infinity)
#define PlusInfinity Infinity

/* Basic operations for type Float */
#define FABS(x)   (x > Zero ? x : -x)
#define FMAX(x,y) (x > y ? x : y)
#define FMIN(x,y) (x > y ? y : x)

/* Setting the IEEE rounding mode bits on a SPARC: */

extern char *ROUND_BUFP;

#define ROUND_DOWN  ieee_flags("set", "direction", "negative", &ROUND_BUFP)
#define ROUND_UP    ieee_flags("set", "direction", "positive", &ROUND_BUFP)
#define ROUND_NEAR  ieee_flags("set", "direction", "nearest",  &ROUND_BUFP)
#define ROUND_ZERO  ieee_flags("set", "direction", "tozero",   &ROUND_BUFP)

/* Procedures: */

void flt_print (FILE *f, Float x);
  /* Prints $x$ on $f$, using "%e" format with full precision. */

/* Error-monitoring arithmetic routines */
/* The routines below perform *zp = x OP y, where OP is the appropriate */
/*   arithmetic operation, using the current rounding mode. */
/* They also add to *errp a quantity d that is an upper bound to */
/*   the error made in the computation of *zp. */
/* They will return *errp = PlusInfinity in case of overflow. */
/* WARNING: they may change the current rounding mode. */

void flt_add(Float x, Float y, FloatP zp, FloatP errp);
  /* Computes $*zp = x + y$.  */

void flt_sub(Float x, Float y, FloatP zp, FloatP errp);
  /* Computes $*zp = x - y$.  */

void flt_mul(Float x, Float y, FloatP zp, FloatP errp);
  /* Computes $*zp = x * y$. */

void flt_inv(Float x, FloatP zp, FloatP errp);
  /* Computes $*zp = 1/x$. */

void flt_sqrt(Float x, FloatP zp, FloatP errp);
  /* Computes $*zp = sqrt(x)$. */

Float flt_random(void);
  /* A random Float value in [0 __ 1) */

#endif
