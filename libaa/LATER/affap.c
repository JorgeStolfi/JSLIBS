/* See affap.h */
/* Last edited on 2017-06-21 05:32:27 by stolfilocal */

#include <affap.h>
#include <affirm.h>
#include <flt.h>
#include <ia.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#define QUICKMUL 1

#define affap_TRACE_ADD      0
#define affap_TRACE_MUL      0
#define affap_TRACE_SQR      0
#define affap_TRACE_NEG      0
#define affap_TRACE_INV      0
#define affap_TRACE_SQRT     0
#define affap_TRACE_EXP      0
#define affap_TRACE_ABS      0
#define affap_TRACE_FIX      0
#define affap_TRACE_THROW    0
#define affap_TRACE_SCALE    0
#define affap_TRACE_MIX      0
#define affap_TRACE_AFFINE   0
#define affap_TRACE_JOIN     0
#define affap_TRACE_CONST    0
#define affap_TRACE_FROM_INT 0

void affap_sqrt
  ( Float xc,
    Float xr,
    FloatP alphap,
    FloatP zetap,
    FloatP gammap,
    FloatP deltap
  )
  /* Computes chebyshev approx to {\sqrt{x}} in {xr}. */
  { Float xlo, xhi, ra, rb, da, db, axa, axb, dlo, dhi;

    #if affap_TRACE_SQRT
      fprintf(stderr, "  enter affap_sqrt:\n");
      fprintf(stderr, "    xc = %15.8e\n", xc);
      fprintf(stderr, "    xr = %15.8e\n", xr);
    #endif

    /* I believe overflow is not a problem except where tested for: */
    (*alphap) = One;
    ROUND_DOWN;
    xlo = xc - xr;
    if (xlo <= Zero) xlo = Zero;
    ra = sqrt(xlo);
    ROUND_UP;
    xhi = xc + xr;
    if(xhi >= PlusInfinity)
      { (*zetap) = One; (*gammap) = Zero; (*deltap) = PlusInfinity; return; }
    rb = sqrt(xhi);
    ROUND_NEAR;
    (*zetap) = ra + rb;
    ROUND_UP;
    axa = ra*ra/(*zetap);
    axb = rb*rb/(*zetap);
    da = axa - ra; da = -da;
    db = axb - rb; db = -db;
    dlo = FMIN(da, db);
    dhi = (*zetap)/Four;
    ROUND_NEAR;
    (*gammap) = (dhi + dlo)/Two;
    ROUND_UP;
    { Float r1 = dhi - (*gammap);
      Float r2 = (*gammap) - dlo;
      (*deltap) = FMAX(r1, r2);
    }
    #if affap_TRACE_SQRT
      fprintf(stderr,
	"    alpha = %15.8e  zeta = %15.8e  gamma  = %15.8e  delta = %15.8e\n",
	*alphap, *zetap, *gammap, *deltap
      );
      fprintf(stderr, "  exit affap_sqrt.\n");
    #endif
  }

void affap_inv(
    Interval xr,
    Float xc,
    FloatP alphap,
    FloatP zetap,
    FloatP gammap,
    FloatP deltap
  )
  /* Computes chebyshev approx to {1/x} in {xr}. */
  {
    Float maxabsalpha, axa, axb, da, db, dlo, dhi;
    int negative;

    #if affap_TRACE_INV
      fprintf(stderr, "  enter affap_inv:\n");
      fprintf(stderr, "    xr.lo = %15.8e\n", xr.lo);
      fprintf(stderr, "    xr.hi = %15.8e\n", xr.hi);
    #endif

    (*zetap) = One;
    /* Assumes {xr} doesn't contain zero: */
    negative = (xr.lo < Zero);
    if (negative)
      { Float t; t = xr.lo; xr.lo = -xr.hi; xr.hi = -t; xc = -xc; }
    /* Compute approximate slope in interval: */
    ROUND_DOWN;
    *alphap = -((One/xr.hi)/xr.lo);
    maxabsalpha = MaxFloat/Two/xc;
    if (-(*alphap) >= maxabsalpha)
      { *alphap = -(maxabsalpha); }
    /* Compute max and min of {1/x - \alpha x} in {xr}: */
    ROUND_NEAR;
    if ( -(*alphap)*xr.hi*xr.hi <= One )
      { /* Difference is monotonically decreasing in {xr}. */
        ROUND_UP;
        axb = (*alphap) * xr.hi;
        ROUND_DOWN;
        axa = (*alphap) * xr.lo;
        dlo = One/xr.hi - axb;
        ROUND_UP;
        dhi = One/xr.lo - axa;
        affirm((dlo <= dhi), "affap_inv: case 1");
      }
    else if ( -(*alphap)*xr.lo*xr.lo >= One )
      { /* Difference is monotonically increasing in {xr}. */
        ROUND_DOWN;
        axb = (*alphap) * xr.hi;
        ROUND_UP;
        axa = (*alphap) * xr.lo;
        dhi = One/xr.hi - axb;
        ROUND_DOWN;
        dlo = One/xr.lo - axa;
        affirm((dlo <= dhi), "affap_inv: case 2");
      }
    else
      { /* Difference may be monotonic or concave in {xr}. */
        ROUND_DOWN;
        dlo = Two * sqrt(-(*alphap));
        axa = (*alphap) * xr.lo;
        axb = (*alphap) * xr.hi;
        ROUND_UP;
        da = One/xr.lo - axa;
        db = One/xr.hi - axb;
        dhi = FMAX(da, db);
        affirm ((dlo <= dhi), "affap_inv: case 3");
      }
    ROUND_NEAR;
    *gammap = (dhi + dlo)/Two;
    ROUND_UP;
    { Float r1 = dhi - (*gammap);
      Float r2 = (*gammap) - dlo;
      *deltap = FMAX(r1, r2);
    }
    if (negative) *gammap = -(*gammap);

    #if affap_TRACE_INV
      fprintf(stderr,
	"    alpha = %15.8e  gamma = %15.8e  delta = %15.8e\n",
	*alphap, *gammap, *deltap
      );
      fprintf(stderr, "  exit affap_inv.\n");
    #endif
  }

void affap_exp(
    Interval xr,
    FloatP alphap,
    FloatP zetap,
    FloatP gammap,
    FloatP deltap
  )
  /* Computes Chebyshev-like approx to {\exp{x}} in {xr}. */
  {
    Float ea, eb, alpha, da, db, dlo, dhi;

    #if affap_TRACE_EXP
      fprintf(stderr, "  enter affap_exp:\n");
      fprintf(stderr, "    xr.lo = %15.8e\n", xr.lo);
      fprintf(stderr, "    xr.hi = %15.8e\n", xr.hi);
    #endif

    (*zetap) = One;
    /* Upper bounds to {\exp(a)} and {\exp(b)}: */
    ROUND_UP;
    eb = exp(xr.hi);
    if (eb >= PlusInfinity)
      { (*alphap) = PlusInfinity;
        (*gammap) = Zero;
        (*deltap) = PlusInfinity;
        return;
      }
    ea = exp(xr.lo);

    /* Select a slope {\alpha}, in or near the range {[\exp(a) _ \exp(b)]} */
    ROUND_NEAR;
    alpha = exp(Half*(xr.lo + xr.hi));
    (*alphap) = alpha;
    
    /* Lower bound for the curve {\exp(x) - \alpha*x}: */
    ROUND_UP;
    dlo = -(alpha*(log(alpha)-1));
    da = ea + (-alpha)*xr.lo;
    db = eb + (-alpha)*xr.hi;
    dhi = FMAX(da, db);
    
    ROUND_NEAR;
    (*gammap) = Half*(dhi + dlo);
    ROUND_UP;
    { Float r1 = dhi - (*gammap);
      Float r2 = (*gammap) - dlo;
      (*deltap) = FMAX(r1, r2);
    }
    #if affap_TRACE_EXP
      fprintf(stderr,
	"    alpha = %15.8e  gamma  = %15.8e  delta = %15.8e\n",
	*alphap, *gammap, *deltap
      );
      fprintf(stderr, "  exit affap_exp.\n");
    #endif
  }

void affap_abs(
    Interval xr,
    FloatP alphap,
    FloatP zetap,
    FloatP gammap,
    FloatP deltap
  )
  /* 
    Computes Chebyshev approximation {\alpha x + \gamma} to {\abs{x}},
    in the zero-straddling interval {xr}. */
  { 
    Float hd, alpha, gamma, delta;

    #if affap_TRACE_ABS
      fprintf(stderr, "  enter affap_abs:\n");
      fprintf(stderr, "    xr.lo = %15.8e\n", xr.lo);
      fprintf(stderr, "    xr.hi = %15.8e\n", xr.hi);
    #endif

    (*zetap) = One;
    ROUND_NEAR; hd = xr.hi/Two - xr.lo/Two; /* Beware of overlow */
    if (hd == Zero)
      { alpha = Zero; }
    else
      { alpha = (xr.hi/Two + xr.lo/Two)/hd; }
    gamma = (One - alpha)*xr.hi; /* Shouldn't overflow */

    /* Estimate approximation error {\delta}: */
    ROUND_UP;
    { Float r1 = (-xr.lo) + alpha*(-xr.lo) - gamma;
      Float r2 = xr.hi + alpha*(-xr.hi) - gamma;
      Float r12 = FMAX(r1, r2);
      delta = FMAX(gamma, r12);
    }
    
    *alphap = alpha;
    *gammap = gamma;
    *deltap = delta;

    #if affap_TRACE_ABS
      fprintf(stderr,
	"    alpha = %15.8e  gamma  = %15.8e  delta = %15.8e\n",
	*alphap, *gammap, *deltap
      );
      fprintf(stderr, "  exit affap_abs.\n");
    #endif
  }
  
