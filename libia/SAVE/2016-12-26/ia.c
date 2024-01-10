/* See ia.h */
/* Last edited on 2016-12-26 17:13:07 by stolfilocal */

#include <affirm.h>
#include <math.h>
#include <stdlib.h>
#include <flt.h>
#include "ia.h"

#define ia_FULL  (ia_Full)
#define ia_ISFULL(x)  (((x).lo <= MinusInfinity) || ((x).hi >= PlusInfinity))
#define ia_NORMFULL(x)  if (ia_ISFULL(x)) (x) = ia_FULL

static Interval ia_Full;

void ia_init(void)
  { flt_init();
    ia_Full.lo = MinusInfinity; ia_Full.hi = PlusInfinity;
  }
  
Interval ia_full(void)
  { return (ia_FULL); }

int ia_is_full(Interval *x)
  { return (ia_ISFULL(*x)); }

void ia_norm_full (Interval *x)
  { ia_NORMFULL(*x); }
  
Interval ia_const(Float x, Float err)
  { Interval z;
    ROUND_DOWN;
    z.lo = x - err;
    ROUND_UP;
    z.hi = x + err;
    ia_NORMFULL(z);
    return z;
  }
    
Interval ia_int_const(int i)
  { Interval z;
    ROUND_DOWN;
    z.lo = flt_from_int(i);
    ROUND_UP;
    z.hi = flt_from_int(i);
    return z;
  }
    
Interval ia_add (Interval x, Interval y)
  { Interval z;
    if (ia_ISFULL(x) || ia_ISFULL(y)) { return ia_FULL; }
    ROUND_DOWN;
    z.lo = x.lo + y.lo;
    ROUND_UP;
    z.hi = x.hi + y.hi;
    ia_NORMFULL(z);
    return z;
  }

Interval ia_sub(Interval x, Interval y)
  { Interval z;
    if (ia_ISFULL(x) || ia_ISFULL(y)) { return ia_FULL; }
    ROUND_DOWN;
    z.lo = x.lo - y.hi;
    ROUND_UP;
    z.hi = x.hi - y.lo;
    ia_NORMFULL(z);
    return z;
  }

Interval ia_neg(Interval x)
  { Interval z;
    if ( ia_ISFULL(x) ) { return ia_FULL; }
    z.lo = -x.hi;
    z.hi = -x.lo;
    return z;
  }

Interval ia_scale (Interval x, Float alpha, Float zeta)
  { if ( ia_ISFULL(x) ) { return ia_FULL; }
    if (alpha == zeta)
      { return(x); }
    else if (alpha == -zeta)
      { return(ia_neg(x)); }
    else if ((zeta == Zero) || (alpha == PlusInfinity))
      { return(ia_FULL); }
    else if ((alpha == Zero) || (zeta == PlusInfinity))
      { return((Interval){Zero, Zero}); }
    else
      {	Interval z;
        double tmp;
        if ((alpha < Zero) == (zeta < Zero))
          { ROUND_DOWN;
	    tmp = (double)alpha * (double)x.lo;
            z.lo = (Float)(tmp/zeta);
	    ROUND_UP;
	    tmp = (double)alpha * (double)x.hi;
	    z.hi = (Float)(tmp/zeta);
          }
        else
          { ROUND_DOWN;
	    tmp = (double)alpha * (double)x.hi;
            z.lo = (Float)(tmp/zeta);
	    ROUND_UP;
	    tmp = (double)alpha * (double)x.lo;
	    z.hi = (Float)(tmp/zeta);
          }
	ia_NORMFULL(z);
	return z;
      }
  }
  
Interval ia_shift (Interval x, Float gamma)
  { Interval z;
    if (ia_ISFULL(x)) { return ia_FULL; }
    ROUND_DOWN;
    z.lo = x.lo + gamma;
    ROUND_UP;
    z.hi = x.hi + gamma;
    ia_NORMFULL(z);
    return z;
  }

Interval ia_mul(Interval x, Interval y)
  { Interval z;
    if (ia_ISFULL(x) || ia_ISFULL(y)) { return ia_FULL; }
    if (x.lo >= Zero)
      {	if (y.lo >= Zero)
	  { ROUND_DOWN;
	    z.lo = x.lo * y.lo;
	    ROUND_UP;
	    z.hi = x.hi * y.hi;
	  }
	else if (y.hi <= Zero)
	  { ROUND_DOWN;
	    z.lo = x.hi * y.lo;
	    ROUND_UP;
	    z.hi = x.lo * y.hi;
	  }
	else
	  { ROUND_DOWN;
	    z.lo = x.hi * y.lo;
	    ROUND_UP;
	    z.hi = x.hi * y.hi;
	  }
      }
    else if (x.hi <= Zero)
      {	if (y.lo >= Zero)
	  { ROUND_DOWN;
	    z.lo = x.lo * y.hi;
	    ROUND_UP;
	    z.hi = x.hi * y.lo;
	  }
	else if (y.hi <= Zero)
	  { ROUND_DOWN;
	    z.lo = x.hi * y.hi;
	    ROUND_UP;
	    z.hi = x.lo * y.lo;
	  }
	else
	  { ROUND_DOWN;
	    z.lo = x.lo * y.hi;
	    ROUND_UP;
	    z.hi = x.lo * y.lo;
	  }
      }
    else
      {	if (y.lo >= Zero)
	  { ROUND_DOWN;
	    z.lo = x.lo * y.hi;
	    ROUND_UP;
	    z.hi = x.hi * y.hi;
	  }
	else if (y.hi <= Zero)
	  { ROUND_DOWN;
	    z.lo = x.hi * y.lo;
	    ROUND_UP;
	    z.hi = x.lo * y.lo;
	  }
	else
	  { ROUND_DOWN;
	    Float z1 = x.lo * y.hi;
            Float z2 = x.hi * y.lo;
            z.lo = FMIN(z1, z2);
	    ROUND_UP;
	    Float z3 = x.lo * y.lo;
            Float z4 = x.hi * y.hi;
            z.hi = FMAX(z3, z4);
	  }
      }
    ia_NORMFULL(z);
    return z;
  }

Interval ia_div (Interval x, Interval y)
  { Interval z = ia_FULL;
    if (ia_ISFULL(x) || ia_ISFULL(y)) { return z; }
    if ((y.lo <= Zero) && (y.hi >= Zero)) { return z; }
    if (x.lo >= Zero)
      {	if (y.lo >= Zero)
	  { ROUND_DOWN;
	    z.lo = x.lo / y.hi;
	    ROUND_UP;
	    z.hi = x.hi / y.lo;
	  }
	else if (y.hi <= Zero)
	  { ROUND_DOWN;
	    z.lo = x.hi / y.hi;
	    ROUND_UP;
	    z.hi = x.lo / y.lo;
	  }
      }
    else if (x.hi <= Zero)
      {	if (y.lo >= Zero)
	  { ROUND_DOWN;
	    z.lo = x.lo / y.lo;
	    ROUND_UP;
	    z.hi = x.hi / y.hi;
	  }
	else if (y.hi <= Zero)
	  { ROUND_DOWN;
	    z.lo = x.hi / y.lo;
	    ROUND_UP;
	    z.hi = x.lo / y.hi;
	  }
      }
    else
      {	if (y.lo >= Zero)
	  { ROUND_DOWN;
	    z.lo = x.lo / y.lo;
	    ROUND_UP;
	    z.hi = x.hi / y.lo;
	  }
	else if (y.hi <= Zero)
	  { ROUND_DOWN;
	    z.lo = x.hi / y.hi;
	    ROUND_UP;
	    z.hi = x.lo / y.hi;
	  }
      }
    ia_NORMFULL(z);
    return z;
  }

Interval ia_sqr(Interval x)
  { Interval z;
    if (ia_ISFULL(x)) { return ia_FULL; }
    if (x.lo >= Zero)
      {	ROUND_DOWN;
	z.lo = x.lo * x.lo;
	ROUND_UP;
	z.hi = x.hi * x.hi;
      }
    else if (x.hi <= Zero)
      {	ROUND_DOWN;
	z.lo = x.hi * x.hi;
	ROUND_UP;
	z.hi = x.lo * x.lo;
      }
    else
      {	z.lo = Zero;
	ROUND_UP;
	Float z1 = x.lo * x.lo;
        Float z2 = x.hi * x.hi;
        z.hi = FMAX(z1, z2);
      }
    ia_NORMFULL(z);
    return z;
  }

Interval ia_inv(Interval x)
  { Interval z = ia_FULL;
    if (ia_ISFULL(x)) { return z; }
    if (x.lo > Zero || x.hi < Zero)
      {	ROUND_DOWN;
	z.lo = One / x.hi;
	ROUND_UP;
	z.hi = One / x.lo;
      }
    else
      { return z; }
    ia_NORMFULL(z);
    return z;
  }

Interval ia_sqrt(Interval x)
  { Interval z;
    if (ia_ISFULL(x)) { return ia_FULL; }
    if (x.hi < Zero)
      {	fatalerror("ia_sqrt: negative interval"); return ia_FULL; }
    else if (x.hi == Zero)
      {	z.hi = Zero; }
    else
      {	ROUND_UP;
	z.hi = (Float)sqrt(x.hi);
      }
    if (x.lo > Zero)
      {	ROUND_DOWN;
	z.lo = (Float)sqrt(x.lo);
      }
    else
      {	z.lo = Zero; }
    return z;
  }

Interval ia_exp(Interval x)
  { Interval z;
    if (ia_ISFULL(x)) { return ia_FULL; }
    if (x.hi == Zero)
      {	z.hi = One; }
    else
      {	ROUND_UP;
	z.hi = (Float)exp(x.hi);
      }
    if (x.lo == Zero)
      {	z.lo = One; }
    else
      {	ROUND_DOWN;
	z.lo = (Float)exp(x.lo);
      }
    return z;
  }

Interval ia_interp (Float x0, Interval y0, Float x1, Interval y1, Float x)
  {
    /* Check for degenerate data: */
    if (x0 == x1)
      { if (x == x0)
          { return ia_join(y0, y1); }
        else
          { return ia_FULL; }
      }
    /* Check for easy cases: */
    if (x == x0)
      { return y0; }
    else if (x == x1)
      { return y1; }
    /* Now we know that {x0 < x1}, and {x} is neither of them. */
    /* Swap if needed so that {x0 < x1}: */
    if (x0 > x1)
      { { Float t = x0; x0 = x1; x1 = t; }
        { Interval t = y0; y0 = y1; y1 = t; }
      }
    /* Check for infinite cases: */
    if (ia_ISFULL(y0) || ia_ISFULL(y1))
      { return ia_FULL; }
    else if ((x0 == MinusInfinity) || (x1 == PlusInfinity))
      { /* Lines are horizontal and must go through both: */
        return ia_meet(y0, y1);
      }
    else if ((x == MinusInfinity) || (x == PlusInfinity))
      { return ia_FULL; }
    /* Interpolate the lines through the proper corners: */
    Interval r; 
    if (x < x0)
      { r.lo = flt_interp_lo(x0, y0.lo, x1, y1.hi, x);
        r.hi = flt_interp_hi(x0, y0.hi, x1, y1.lo, x);
      }
    else if (x > x1)
      { r.lo = flt_interp_lo(x0, y0.hi, x1, y1.lo, x);
        r.hi = flt_interp_hi(x0, y0.lo, x1, y1.hi, x);
      }
    else 
      { r.lo = flt_interp_lo(x0, y0.lo, x1, y1.lo, x);
        r.hi = flt_interp_hi(x0, y0.hi, x1, y1.hi, x);
      }
    return r;
  }

Interval ia_affine(Interval x, Float alpha, Float zeta, Float gamma, Float delta)
  {
    delta = FABS(delta);
    if (
      ia_ISFULL(x) ||
      (FABS(alpha) >= PlusInfinity) ||
      (FABS(zeta) == Zero) ||
      (FABS(gamma) >= PlusInfinity) ||
      (delta >= PlusInfinity)
    ) 
      { return (ia_FULL); }
    else
      { Interval z;

        /* Ensure that alpha and zeta are non-negative: */

        if (zeta < Zero)
          { alpha = -alpha; zeta = -zeta; }
        if (alpha < Zero) 
          { Float t = x.lo; x.lo = -x.hi; x.hi = -t; alpha = -alpha; }
          
        /* Compute result: */
        
        if ((alpha == Zero) || (FABS(zeta) >= PlusInfinity))
	  { ROUND_DOWN;
	    z.lo = gamma - delta;
	    ROUND_UP;
	    z.hi = gamma + delta;
	  }
	else if (alpha == zeta)
	  { ROUND_DOWN;
	    z.lo = x.lo + gamma - delta;
	    ROUND_UP;
	    z.hi = x.hi + gamma + delta;
	  }
	else 
	  { double tmp;
	    ROUND_DOWN;
	    tmp = (double)alpha * (double)x.lo;
	    z.lo = ((Float)tmp / zeta) + gamma - delta;
	    ROUND_UP;
	    tmp = (double)alpha * (double)x.hi;
	    z.hi = ((Float)tmp / zeta) + gamma + delta;
	  }
	ia_NORMFULL(z);
	return z;
      }
  }

Interval ia_affine_2(
    Interval x, Float alpha,
    Interval y, Float beta, 
    Float zeta, 
    Float gamma, Float delta
  )
  {
    delta = FABS(delta);
    if (
      ia_ISFULL(x) ||
      (FABS(alpha) >= PlusInfinity) ||
      (FABS(beta) >= PlusInfinity) ||
      (FABS(zeta) == Zero) ||
      (FABS(gamma) >= PlusInfinity) ||
      (delta >= PlusInfinity)
    ) 
      { return (ia_FULL); }
    else
      { Interval z;
	
        /* Ensure that alpha, beta, and gamma are non-negative: */

        if (zeta < Zero)
          { alpha = -alpha; beta = -beta; zeta = -zeta; }
        if (alpha < Zero) 
          { Float t = x.lo; x.lo = -x.hi; x.hi = -t; alpha = -alpha; }
        if (beta < Zero )
          { Float t = y.lo; y.lo = -y.hi; y.hi = -t; beta = -beta; }
          
        /* Compute result: */

        if (FABS(zeta) >= PlusInfinity)
	  { ROUND_DOWN;
	    z.lo = gamma - delta;
	    ROUND_UP;
	    z.hi = gamma + delta;
	  }
        else if (alpha == Zero)
          { return(ia_affine(y, beta, zeta, gamma, delta)); }
        else if (beta == Zero)
          { return(ia_affine(x, alpha, zeta, gamma, delta)); }
        else
	  { double tmp;
	    ROUND_DOWN;
	    tmp = (double)alpha * (double)x.lo + (double)beta * (double)y.lo;
	    z.lo = ((Float)tmp / zeta) + gamma - delta;
	    ROUND_UP;
	    tmp = (double)alpha * (double)x.hi + (double)beta * (double)y.hi;
	    z.hi = ((Float)tmp / zeta) + gamma + delta;
	  }
	ia_NORMFULL(z);
	return z;
      }
  }

Interval ia_abs   (Interval x)
  { if (x.lo >= Zero)
      { return (x); }
    else if (x.hi <= Zero) 
      { return ((Interval){-x.hi, -x.lo}); }
    else 
      { return ((Interval){Zero, FMAX((-x.lo), x.hi)}); }
  }
  
Interval ia_max   (Interval x, Interval y)
  { return ((Interval){FMAX(x.lo, y.lo), FMAX(x.hi, y.hi)}); }
  
Interval ia_min   (Interval x, Interval y)
  { return ((Interval){FMIN(x.lo, y.lo), FMIN(x.hi, y.hi)}); }

Interval ia_throw (void)
  { int coins;
    coins = rand();
    if ((coins&255) == 0)
      return (ia_full());
    else if ((coins&63) == 0)
      return ((Interval){Zero, Zero});
    else
      { Interval z;
	Float ctr = flt_random_mag(0, 25) * flt_random();
	if ((coins & 32) == 0) ctr = -ctr;
        if ((coins & 3) == 0)
          { z.lo = FMIN(ctr, Zero);
            z.hi = FMAX(ctr, Zero);
          }
        else
          { Float rad = flt_random_mag(0, 25) * flt_random();
            ROUND_DOWN; z.lo = ctr - rad;
	    ROUND_UP; z.hi = ctr + rad;
          }
	if ((z.lo <= MinusInfinity) || (z.hi >= PlusInfinity)) { z = ia_full(); }
        return(z);
      }
  }
    
void ia_print_bound(FILE *f, Float v, int which, int full)
  { int i;
    if (full)
      { putc(' ', f);
        for (i = 1; i < flt_FMT_WIDTH; i++) putc('*', f);
      }
    else
      { if (which == 0) { ROUND_DOWN; } else { ROUND_UP; }
        fprintf(f, flt_FMT_SPEC, v);
      }
  }

void ia_print (FILE *f, Interval x)
  { int full = ia_ISFULL(x);
    putc('[', f);
    ia_print_bound(f, x.lo, 0, full);
    fputs(" __ ", f);
    ia_print_bound(f, x.hi, 1, full);
    putc(']', f);
  }

Float ia_mid (Interval x)
  { if (ia_ISFULL(x))
      { return 0; }
    else if (x.lo == x.hi)
      { return x.lo; }
    else
      { Float m;
        ROUND_NEAR;
        m = (Float)((x.lo * Half) + (x.hi * Half));
        affirm((m >= x.lo) && (m <= x.hi), "mean failed");
        return m;
      }
  }

Float ia_rad (Interval x)
  { if (ia_ISFULL(x))
      { return PlusInfinity; }
    else if (x.lo == x.hi)
      { return Zero; }
    else
      { Float m, rlo, rhi;
        ROUND_NEAR;
        m = (Float)((x.lo * Half) + (x.hi * Half));
        rlo = (Float)(m - x.lo);
        rhi = (Float)(x.hi - m);
        affirm((rlo >= Zero) && (rhi >= Zero), "radius failed");
        return (rlo > rhi ? rlo: rhi);
      }
  }

Interval ia_meet(Interval x, Interval y)
  { Interval z;
    z.lo = FMAX(x.lo, y.lo);
    z.hi = FMIN(x.hi, y.hi);
    return z;
  }

Interval ia_join (Interval x, Interval y)
  { Interval z;
    z.lo = FMIN(x.lo, y.lo);
    z.hi = FMAX(x.hi, y.hi);
    return z;
  }
   
