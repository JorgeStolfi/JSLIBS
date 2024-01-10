/* first-order interval arithmetic routines */
/* Jorge Stolfi 13 jan 1993                 */

#include "interval.h"
#include "foifloat.h"
#include "foimisc.h"
#include "iomisc.h"
#include <math.h>

Interval iv_full (void)
  { Interval z;
    z.lo = MinusInfinity; z.hi = PlusInfinity;
    return(z);
  }
  
Interval iv_const(Float x, Float err)
  {
    Interval z;
    ROUND_DOWN;
    z.lo = x - err;
    ROUND_UP;
    z.hi = x + err;
    IV_NORMFULL(z);
    return (z);
  }
    
Interval iv_add (Interval x, Interval y)
  {
    Interval z;
    if (IV_ISFULL(x) || IV_ISFULL(y)) return (IV_FULL);
    ROUND_DOWN;
    z.lo = x.lo + y.lo;
    ROUND_UP;
    z.hi = x.hi + y.hi;
    IV_NORMFULL(z);
    return (z);
  }

Interval iv_neg(Interval x)
  {
    Interval z;
    if ( IV_ISFULL(x) ) return (IV_FULL);
    z.lo = -x.hi;
    z.hi = -x.lo;
    return (z);
  }

Interval iv_sub(Interval x, Interval y)
  {
    Interval z;
    if (IV_ISFULL(x) || IV_ISFULL(y)) return (IV_FULL);
    ROUND_DOWN;
    z.lo = x.lo - y.hi;
    ROUND_UP;
    z.hi = x.hi - y.lo;
    IV_NORMFULL(z);
    return (z);
  }

Interval iv_mul(Interval x, Interval y)
  {
    Interval z;
    if (IV_ISFULL(x) || IV_ISFULL(y)) return (IV_FULL);
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
    IV_NORMFULL(z);
    return (z);
  }

Interval iv_sqr(Interval x)
  {
    Interval z;
    if (IV_ISFULL(x)) return (IV_FULL);
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
    IV_NORMFULL(z);
    return (z);
  }

Interval iv_inv(Interval x)
  {
    Interval z;
    if (IV_ISFULL(x)) return (IV_FULL);
    if (x.lo > Zero || x.hi < Zero)
      {
	ROUND_DOWN;
	z.lo = One / x.hi;
	ROUND_UP;
	z.hi = One / x.lo;
      }
    else
      return (IV_FULL);
    IV_NORMFULL(z);
    return (z);
  }

Interval iv_sqrt(Interval x)
  {
    Interval z;
    if (IV_ISFULL(x)) return (IV_FULL);
    if (x.hi < Zero)
      {
	error("iv_sqrt: negative interval");
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

void iv_print (FILE *f, Interval x)
  {
    int i;
    putc('[', f);
    if (IV_ISFULL(x))
      {
        for (i=0;i<F_FMT_WIDTH;i++) putc('*', f);
        fputs(" _ ", f);
        for (i=0;i<F_FMT_WIDTH;i++) putc('*', f);
      }
    else
      {
        ROUND_DOWN;
        fprintf(f, F_FMT_SPEC, x.lo);
        fputs(" _ ", f);
        ROUND_UP;
        fprintf(f, F_FMT_SPEC, x.hi);
      }
    putc(']', f);
  }

Interval iv_meet(Interval x, Interval y)
  { Interval z;
    z.lo = FMAX(x.lo, y.lo);
    z.hi = FMIN(x.hi, y.hi);
    if (z.lo > z.hi)
      { iv_print(stderr, x); putc('\n', stderr);
        iv_print(stderr, y); putc('\n', stderr);
        error ("iv_meet: disjoint intervals");
      }
    return (z);
  }

Interval iv_affine(
    Interval x,
    Float alpha,
    Float beta,
    Float gamma
  )
  {
    Interval z;
    gamma = FABS(gamma);
    if (IV_ISFULL(x)
    || (FABS(alpha) >= PlusInfinity)
    || (FABS(beta) >= PlusInfinity)
    || (gamma >= PlusInfinity)
    ) return (IV_FULL);
    if (alpha == Zero)
      { ROUND_DOWN;
        z.lo = beta - gamma;
        ROUND_UP;
        z.hi = beta + gamma;
      }
    else if (alpha < Zero)
      { ROUND_DOWN;
        z.lo = alpha * x.hi + beta - gamma;
        ROUND_UP;
        z.hi = alpha * x.lo + beta + gamma;
      }
    else
      { ROUND_DOWN;
        z.lo = alpha * x.lo + beta - gamma;
        ROUND_UP;
        z.hi = alpha * x.hi + beta + gamma;
      }
    IV_NORMFULL(z);
    return (z);
  }

