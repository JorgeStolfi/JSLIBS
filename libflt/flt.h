/* Operations on floating-point values */
/* Last edited on 2016-04-01 02:08:49 by stolfilocal */

#ifndef flt_H
#define flt_H

#define _GNU_SOURCE
#include <stdio.h>
/* #include <sys/ieeefp.h> */

/* The floating-point type used by the operations below.: */
typedef
  float Float;

typedef
  float *FloatP;

/* Default {printf} parameters for type Float: */

#define flt_FMT_DIGITS 8
#define flt_FMT_WIDTH (flt_FMT_DIGITS + 6)
#define flt_FMT_DECIMALS (flt_FMT_DIGITS - 1)
#define flt_FMT_SPEC "%14.7e"

/* Generator for a random Float number in [0 _ 1): */

#define flt_RANDOM frandom

/* Common constants for type Float: */

#define One     (1.0f)
#define Zero    (0.0f)
#define Half    (0.5f)
#define Quarter (0.25f)
#define Two     (2.0f)
#define Three   (3.0f)
#define Four    (4.0f)
#define Eight   (8.0f)
#define Sixteen (16.0f)

extern Float MaxFloat;
extern Float Infinity;

#define MinusInfinity (-Infinity)
#define PlusInfinity Infinity

/* Basic operations for type {Float} */
#define FABS(x)   ((x) > Zero ? (x) : -(x))
#define FMAX(x,y) ((x) > (y) ? (x) : (y))
#define FMIN(x,y) ((x) > (y) ? (y) : (x))

/* Procedures: */

void flt_init (void);
  /* Initializes the constants ({Infinity}, {MaxFloat}, {MinusInfinity}, etc.).
    This routine MUST be called at least once before any other routine below. */

void flt_print (FILE *f, Float x);
  /* Prints {x} on {f}, using "%e" format with full precision. */

/* 
  Error-monitoring arithmetic routines
  
  The routines below perform {*zp = x OP y}, where {OP} is the appropriate
  arithmetic operation, using the current rounding mode.
  
  They also add to {*errp} (with upwards rounding) a quantity that is
  an upper bound to the error made in the computation of {*zp}.
  
  They will return {*errp = PlusInfinity} in case of overflow.
  WARNING: they may change the current rounding mode. */

void flt_add(Float x, Float y, FloatP zp, FloatP errp);
  /* Computes {*zp = x + y}.  */

void flt_sub(Float x, Float y, FloatP zp, FloatP errp);
  /* Computes {*zp = x - y}.  */

void flt_mul(Float x, Float y, FloatP zp, FloatP errp);
  /* Computes {*zp = x * y}. */

void flt_div(Float x, Float y, FloatP zp, FloatP errp);
  /* Computes {*zp = x / y}. */

void flt_inv(Float x, FloatP zp, FloatP errp);
  /* Computes {*zp = 1/x}. */

void flt_sqrt(Float x, FloatP zp, FloatP errp);
  /* Computes {*zp = sqrt(x)}. */

void flt_exp(Float x, FloatP zp, FloatP errp);
  /* Computes {*zp = exp(x)}. */

void flt_scale (
    Float x, Float alpha, 
    Float zeta,
    FloatP zp,
    FloatP errp
  );
  /* Computes {*zp = \alpha x / zeta}. */

void flt_mix (
    Float x, Float alpha, 
    Float y, Float beta, 
    Float zeta,
    FloatP zp,
    FloatP errp
  );
  /* Computes {*zp = ( \alpha x + \beta y ) / \zeta}. */

Float flt_random(void);
  /* A random Float value in [0 __ 1).
     The client must have called {srandom(<seed>)}. */

Float flt_random_mag(int avg, int dev);
  /* Two raised to a random exponent, uniform between 
    {avg-dev} and {avg+dev}.  (The result may be infinity). 
    The client must have called {srandom(<seed>)}. */

Float flt_from_int(int i);
  /* Converts the integer {i} to a {Float} value, according to 
    the current rounding mode. */

/* Setting the IEEE rounding mode bits on a SPARC: */

#define ROUND_DOWN  flt_round_down()
#define ROUND_UP    flt_round_up()
#define ROUND_NEAR  flt_round_near()
#define ROUND_ZERO  flt_round_zero()

void flt_round_down (void);
void flt_round_up   (void);
void flt_round_near (void);
void flt_round_zero (void);
  /* These routines set the IEEE FP rounding direction. */
  /* They should be a lot faster than "ieee_flags". */

int flt_get_fsr(void);
  /* Returns the FP status register. (For debugging) */
  
void flt_single_prec (void);
void flt_double_prec (void);
void flt_extended_prec (void);
  /* Changes the result of arithmetic ops to SINGLE, DOUBLE, EXTENDED */

Float flt_mix_lo(Float xa, Float wa, Float xb, Float wb);
Float flt_mix_hi(Float xa, Float wa, Float xb, Float wb);
  /* Return lower and upper bounds for the weighted average of {xa}
    and {xb} with weights {wa} and {wb} --- i.e. for the formula
    {(xa*wa + xb*wb)/(wa+wb)}. Requires {wa>=0} and {wb>=0}. */

Float flt_interp_lo(Float xa, Float ya, Float xb, Float yb, Float x);
Float flt_interp_hi(Float xa, Float ya, Float xb, Float yb, Float x);
  /* These procedures compute upper and lower bounds, respectively,
    for the interpolation formula {((xb-x)*ya + (x-xa)*yb)/(xb-xa)}
    --- i.e. the Y value at abscissa {x} of the straight line through
    the points {(xa,ya)} and {(xb,yb)}. Both require {xa<=x<=xb}. */
    
#endif
