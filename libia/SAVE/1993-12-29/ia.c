/* See ia.h */

#include "ia.h"
#include <flt.h>
#include <js.h>
#include <ioprotos.h>
#include <math.h>

#define IA_FULL  (ia_Full)
#define IA_ISFULL(x)  (((x).lo <= MinusInfinity) || ((x).hi >= PlusInfinity))
#define IA_NORMFULL(x)  if (IA_ISFULL(x)) (x) = IA_FULL

static Interval ia_Full;

void ia_init(void)
  { flt_init();
    ia_Full.lo = MinusInfinity; ia_Full.hi = PlusInfinity;
  }
  
Interval ia_full(void)
  { return (IA_FULL); }

int ia_is_full(Interval *x)
  { return (IA_ISFULL(*x)); }

void ia_norm_full (Interval *x)
  { IA_NORMFULL(*x); }

  
Interval ia_const(Float x, Float err)
  {
    Interval z;
    ROUND_DOWN;
    z.lo = x - err;
    ROUND_UP;
    z.hi = x + err;
    IA_NORMFULL(z);
    return (z);
  }
    
Interval ia_int_const(int i)
  {
    Interval z;
    ROUND_DOWN;
    z.lo = (Float) i;
    ROUND_UP;
    z.hi = (Float) i;
    return (z);
  }
    
Interval ia_add (Interval x, Interval y)
  {
    Interval z;
    if (IA_ISFULL(x) || IA_ISFULL(y)) return (IA_FULL);
    ROUND_DOWN;
    z.lo = x.lo + y.lo;
    ROUND_UP;
    z.hi = x.hi + y.hi;
    IA_NORMFULL(z);
    return (z);
  }

Interval ia_sub(Interval x, Interval y)
  {
    Interval z;
    if (IA_ISFULL(x) || IA_ISFULL(y)) return (IA_FULL);
    ROUND_DOWN;
    z.lo = x.lo - y.hi;
    ROUND_UP;
    z.hi = x.hi - y.lo;
    IA_NORMFULL(z);
    return (z);
  }

Interval ia_neg(Interval x)
  {
    Interval z;
    if ( IA_ISFULL(x) ) return (IA_FULL);
    z.lo = -x.hi;
    z.hi = -x.lo;
    return (z);
  }

Interval ia_scale (Interval x, Float alpha, Float zeta)
  {
    if ( IA_ISFULL(x) ) return (IA_FULL);
    if (alpha == zeta)
      { return(x); }
    else if (alpha == -zeta)
      { return(ia_neg(x)); }
    else if ((zeta == Zero) || (alpha == PlusInfinity))
      { return(IA_FULL); }
    else if ((alpha == Zero) || (zeta == PlusInfinity))
      { return((Interval){Zero, Zero}); }
    else
      {
	Interval z;
        double tmp;
        if ((alpha < Zero) == (zeta < Zero))
          { ROUND_DOWN;
	    tmp = (double)alpha * (double)x.lo;
            z.lo = tmp/zeta ;
	    ROUND_UP;
	    tmp = (double)alpha * (double)x.hi;
	    z.hi = tmp/zeta;
          }
        else
          { ROUND_DOWN;
	    tmp = (double)alpha * (double)x.hi;
            z.lo = tmp/zeta ;
	    ROUND_UP;
	    tmp = (double)alpha * (double)x.lo;
	    z.hi = tmp/zeta;
          }
	IA_NORMFULL(z);
	return (z);
      }
  }
  
Interval ia_shift (Interval x, Float gamma)
  {
    Interval z;
    if (IA_ISFULL(x)) return (IA_FULL);
    ROUND_DOWN;
    z.lo = x.lo + gamma;
    ROUND_UP;
    z.hi = x.hi + gamma;
    IA_NORMFULL(z);
    return (z);
  }

Interval ia_mul(Interval x, Interval y)
  {
    Interval z;
    if (IA_ISFULL(x) || IA_ISFULL(y)) return (IA_FULL);
    if (x.lo >= Zero)
      {
	if (y.lo >= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.lo * y.lo;
	    ROUND_UP;
	    z.hi = x.hi * y.hi;
	  }
	else if (y.hi <= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.hi * y.lo;
	    ROUND_UP;
	    z.hi = x.lo * y.hi;
	  }
	else
	  {
	    ROUND_DOWN;
	    z.lo = x.hi * y.lo;
	    ROUND_UP;
	    z.hi = x.hi * y.hi;
	  }
      }
    else if (x.hi <= Zero)
      {
	if (y.lo >= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.lo * y.hi;
	    ROUND_UP;
	    z.hi = x.hi * y.lo;
	  }
	else if (y.hi <= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.hi * y.hi;
	    ROUND_UP;
	    z.hi = x.lo * y.lo;
	  }
	else
	  {
	    ROUND_DOWN;
	    z.lo = x.lo * y.hi;
	    ROUND_UP;
	    z.hi = x.lo * y.lo;
	  }
      }
    else
      {
	if (y.lo >= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.lo * y.hi;
	    ROUND_UP;
	    z.hi = x.hi * y.hi;
	  }
	else if (y.hi <= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.hi * y.lo;
	    ROUND_UP;
	    z.hi = x.lo * y.lo;
	  }
	else
	  {
	    ROUND_DOWN;
	    z.lo = FMIN(x.lo * y.hi, x.hi * y.lo);
	    ROUND_UP;
	    z.hi = FMAX(x.lo * y.lo, x.hi * y.hi);
	  }
      }
    IA_NORMFULL(z);
    return (z);
  }

Interval ia_div (Interval x, Interval y)
  {
    Interval z;
    if (IA_ISFULL(x) || IA_ISFULL(y)) return (IA_FULL);
    if ((y.lo <= Zero) && (y.hi >= Zero)) return (IA_FULL);
    if (x.lo >= Zero)
      {
	if (y.lo >= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.lo / y.hi;
	    ROUND_UP;
	    z.hi = x.hi / y.lo;
	  }
	else if (y.hi <= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.hi / y.hi;
	    ROUND_UP;
	    z.hi = x.lo / y.lo;
	  }
      }
    else if (x.hi <= Zero)
      {
	if (y.lo >= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.lo / y.lo;
	    ROUND_UP;
	    z.hi = x.hi / y.hi;
	  }
	else if (y.hi <= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.hi / y.lo;
	    ROUND_UP;
	    z.hi = x.lo / y.hi;
	  }
      }
    else
      {
	if (y.lo >= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.lo / y.lo;
	    ROUND_UP;
	    z.hi = x.hi / y.lo;
	  }
	else if (y.hi <= Zero)
	  {
	    ROUND_DOWN;
	    z.lo = x.hi / y.hi;
	    ROUND_UP;
	    z.hi = x.lo / y.hi;
	  }
      }
    IA_NORMFULL(z);
    return (z);
  }

Interval ia_sqr(Interval x)
  {
    Interval z;
    if (IA_ISFULL(x)) return (IA_FULL);
    if (x.lo >= Zero)
      {
	ROUND_DOWN;
	z.lo = x.lo * x.lo;
	ROUND_UP;
	z.hi = x.hi * x.hi;
      }
    else if (x.hi <= Zero)
      {
	ROUND_DOWN;
	z.lo = x.hi * x.hi;
	ROUND_UP;
	z.hi = x.lo * x.lo;
      }
    else
      {
	z.lo = Zero;
	ROUND_UP;
	z.hi = FMAX(x.lo * x.lo, x.hi * x.hi);
      }
    IA_NORMFULL(z);
    return (z);
  }

Interval ia_inv(Interval x)
  {
    Interval z;
    if (IA_ISFULL(x)) return (IA_FULL);
    if (x.lo > Zero || x.hi < Zero)
      {
	ROUND_DOWN;
	z.lo = One / x.hi;
	ROUND_UP;
	z.hi = One / x.lo;
      }
    else
      return (IA_FULL);
    IA_NORMFULL(z);
    return (z);
  }

Interval ia_sqrt(Interval x)
  {
    Interval z;
    if (IA_ISFULL(x)) return (IA_FULL);
    if (x.hi < Zero)
      {
	error("ia_sqrt: negative interval");
      }
    else if (x.hi == Zero)
      {
	z.hi = Zero;
      }
    else
      {
	ROUND_UP;
	z.hi = sqrt(x.hi);
      }

    if (x.lo > Zero)
      {
	ROUND_DOWN;
	z.lo = sqrt(x.lo);
      }
    else
      {
	z.lo = Zero;
      }

    return (z);
  }

Interval ia_affine (
    Interval x, Float alpha, 
    Float zeta,
    Float gamma, Float delta
  )
  {
    delta = FABS(delta);
    if (
      IA_ISFULL(x) ||
      (FABS(alpha) >= PlusInfinity) ||
      (FABS(zeta) == Zero) ||
      (FABS(gamma) >= PlusInfinity) ||
      (delta >= PlusInfinity)
    ) 
      { return (IA_FULL); }
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
	IA_NORMFULL(z);
	return (z);
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
      IA_ISFULL(x) ||
      (FABS(alpha) >= PlusInfinity) ||
      (FABS(beta) >= PlusInfinity) ||
      (FABS(zeta) == Zero) ||
      (FABS(gamma) >= PlusInfinity) ||
      (delta >= PlusInfinity)
    ) 
      { return (IA_FULL); }
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
	IA_NORMFULL(z);
	return (z);
      }
  }

Interval ia_abs   (Interval x)
  { 
    if (x.lo >= Zero)
      { return (x); }
    else if (x.hi <= Zero) 
      { return ((Interval){-x.hi, -x.lo}); }
    else 
      { return ((Interval){Zero, FMAX(-x.lo, x.hi)}); }
  }
  
Interval ia_max   (Interval x, Interval y)
  { 
    return ((Interval){FMAX(x.lo, y.lo), FMAX(x.hi, y.hi)});
  }
  
Interval ia_min   (Interval x, Interval y)
  {
    return ((Interval){FMIN(x.lo, y.lo), FMIN(x.hi, y.hi)});
  }

void ia_print (FILE *f, Interval x)
  {
    int i;
    putc('[', f);
    if (IA_ISFULL(x))
      {
        putc(' ', f);
        for (i=1;i<FLT_FMT_WIDTH;i++) putc('*', f);
        fputs(" __ ", f);
        putc(' ', f);
        for (i=1;i<FLT_FMT_WIDTH;i++) putc('*', f);
      }
    else
      {
        ROUND_DOWN;
        fprintf(f, FLT_FMT_SPEC, x.lo);
        fputs(" __ ", f);
        ROUND_UP;
        fprintf(f, FLT_FMT_SPEC, x.hi);
      }
    putc(']', f);
  }

Interval ia_meet(Interval x, Interval y)
  { Interval z;
    z.lo = FMAX(x.lo, y.lo);
    z.hi = FMIN(x.hi, y.hi);
    if (z.lo > z.hi)
      { ia_print(stderr, x); putc('\n', stderr);
        ia_print(stderr, y); putc('\n', stderr);
        error ("ia_meet: disjoint intervals");
      }
    return (z);
  }

Interval ia_join (Interval x, Interval y)
  { Interval z;
    z.lo = FMIN(x.lo, y.lo);
    z.hi = FMAX(x.hi, y.hi);
    return (z);
  }
   
