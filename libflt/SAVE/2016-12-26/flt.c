/* See flt.h */
/* Last edited on 2016-04-01 02:44:16 by stolfilocal */

#define _GNU_SOURCE
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <affirm.h>
#include <jsrandom.h>

#include "flt.h"

char *ROUND_BUFP = NULL;

Float MaxFloat;
Float Infinity;

void flt_init(void)
  { Infinity = INFINITY;
    MaxFloat = FLT_MAX;
    ROUND_NEAR; 
    /* fprintf(stderr, "flt_init: fp_precision = %8x  fsr = %8x\n", fp_precision, flt_get_fsr()); */
  }

void flt_print (FILE *f, Float x)
  { fprintf(f, flt_FMT_SPEC, x); }
  
/* There should be a simpler and faster way to do these: */

#define UPDATE_ERR                                                 \
    if (zlo != zhi)                                                \
      {                                                            \
	if ((zhi >= PlusInfinity) || (zlo <= MinusInfinity))       \
	  { *zp = PlusInfinity; *errp = PlusInfinity; return; }    \
	ROUND_UP;                                                  \
	{ Float d1 = zhi - (*zp);                                  \
	  Float d2 = (*zp) - zlo;                                  \
	  Float del = FMAX(d1, d2); /* Shouldn't overflow. */      \
	  *errp = *errp + del;  /* May overflow */                 \
	}                                                          \
      }
  
void flt_add(Float x, Float y, FloatP zp, FloatP errp)
  { Float zhi, zlo;
    *zp = x + y;         /* May overflow. */
    ROUND_DOWN;
    zlo = x + y;
    ROUND_UP;
    zhi = x + y;
    UPDATE_ERR
  }

void flt_sub(Float x, Float y, FloatP zp, FloatP errp)
  { Float zhi, zlo;
    *zp = x - y;         /* May overflow. */
    ROUND_DOWN;
    zlo = x - y;
    ROUND_UP;
    zhi = x - y;
    UPDATE_ERR
  }

void flt_mul(Float x, Float y, FloatP zp, FloatP errp)
  { Float zhi, zlo;
    *zp = x * y;
    ROUND_DOWN;
    zlo = x * y;
    ROUND_UP;
    zhi = x * y;
    UPDATE_ERR
  }

void flt_div(Float x, Float y, FloatP zp, FloatP errp)
  { Float zhi, zlo;
    *zp = x / y;         /* May overflow. */
    ROUND_DOWN;
    zlo = x / y;
    ROUND_UP;
    zhi = x / y;
    UPDATE_ERR
  }

void flt_inv(Float x, FloatP zp, FloatP errp)
  { Float zhi, zlo;
    *zp = One/x;
    ROUND_DOWN;
    zlo = One/x;
    ROUND_UP;
    zhi = One/x;
    UPDATE_ERR
  }

void flt_sqrt(Float x, FloatP zp, FloatP errp)
  { Float zhi, zlo;
    affirm(x >= Zero, "flt_sqrt: argument is negative");
    *zp = (Float)sqrt(x);
    ROUND_DOWN;
    zlo = (Float)sqrt(x);
    ROUND_UP;
    zhi = (Float)sqrt(x);
    UPDATE_ERR;
  }

void flt_exp(Float x, FloatP zp, FloatP errp)
  { Float zhi, zlo;
    affirm(x >= Zero, "flt_sqrt: argument is negative");
    *zp = (Float)exp(x);
    ROUND_DOWN;
    zlo = (Float)exp(x);
    ROUND_UP;
    zhi = (Float)exp(x);
    UPDATE_ERR;
  }

void flt_scale (
    Float x, Float alpha, 
    Float zeta,
    FloatP zp,
    FloatP errp
  )
  { /* Make sure that {zeta} is not negative: */
    if (zeta < 0) { alpha = -alpha; zeta = -zeta; }
    /* Compute the exact product {x*alpha}: */
    double xa = ((double) x) * ((double) alpha); 
    /* With IEEE arithmetic, {xa} should be exact        */
    /*   and without overflow, whatever the current rounding mode. */

    Float zhi, zlo;
    if (zeta == One)
      { *zp = (Float) xa;
	ROUND_DOWN;
	zlo = (Float) xa;
	ROUND_UP;
	zhi = (Float) xa;
      }
    else
      { *zp = ((Float) xa) / zeta;
        ROUND_DOWN;
        zlo = ((Float) xa) / zeta;
        ROUND_UP;
        zhi = ((Float) xa) / zeta;
      }
    UPDATE_ERR
  }

void flt_mix (
    Float x, Float alpha, 
    Float y, Float beta, 
    Float zeta,
    FloatP zp,
    FloatP errp
  )
  { /* Make sure that {zeta} is not negative: */
    if (zeta < 0) { alpha = -alpha; beta = -beta; zeta = -zeta; }
    /* Compute the exact mix {x*alpha + y*beta}: */
    double xa = ((double) x) * ((double) alpha); 
    double yb = ((double) y) * ((double) beta);  
    /* With IEEE arithmetic,  {xa} and {yb} should be exact 
      and without overflow, whatever the current rounding mode. */
    Float zhi, zlo;
    if (zeta == One)
      {  *zp = (Float)(xa + yb);
	ROUND_DOWN;
	zlo = (Float)(xa + yb);
	ROUND_UP;
	zhi = (Float)(xa + yb);
      }
    else
      { *zp = ((Float)(xa + yb)) / zeta;
	ROUND_DOWN;
	zlo = ((Float)(xa + yb)) / zeta;
	ROUND_UP;
	zhi = ((Float)(xa + yb)) / zeta;
      }
    UPDATE_ERR
  }
            
Float flt_random(void)
  {
    return flt_RANDOM();
  }

Float flt_random_mag(int avg, int dev)
  {
#   define TwoToSixteen (65536.0f)
    Float t = One;
    int exp;

    if (dev == 0)
      { exp = avg; }
    else
      { /* Compute exponent "exp": */
	int w = dev + dev + 1;
        unsigned mask = 1;
	/* find smallest all-ones mask that is no less than 2*dev + 1: */
	while ((mask & w) != w) mask = ((mask << 1) | 1);
	/* Generate random integer in range [0..2*dev]: */
	do { exp = (rand() & mask); } while (exp >= w);
	/* Compute exponent: */
	exp = avg + (exp - dev);
      }
    
    /* Compute power: */
    { int j = 0;
      while (j + 16 <= exp) { t *= TwoToSixteen; j += 16; }
      while (j < exp) { t *= Two; j++; }
      while (j - 16 >= exp) { t /= TwoToSixteen; j -= 16; }
      while (j > exp) { t /= Two; j--; }
    }
    
    return(t);
#   undef TwoToSixteen
  }
  
Float flt_from_int(int i)
  {
    double di = (double)i; /* Should be exact! */
    affirm((int)di == i, "bug in double/int conversion");
    Float res = (Float) di;
    return res;
  }

Float flt_interp_lo(Float xa, Float ya, Float xb, Float yb, Float x)
  { demand (xa < xb, "data points out of order"); 
    if (x == xa) 
      { return(ya); }
    else if (x == xb)
      { return(yb); }
    else 
      { ROUND_DOWN;
        /* Choose a reference point {(xref,yref)}: */
        Float xref, yref;
        if (x <= xa)
          { xref = xa; yref = ya; }
        else if (x >= xb)
          { xref = xb; yref = yb; }
        else if (FABS(x - xa) < FABS(x - xb))
          { xref = xa; yref = ya; }
        else 
          { xref = xb; yref = yb; }

        /* Compute the slope and displacement, rounded appropriately: */
        Float slope;
        if (x <= xref) 
          { /* We must round the slope up: */
            slope = -((ya - yb)/(xb - xa));
           }
        else
          { /* We must round the slope down: */
            slope = (yb - ya)/(-(xa - xb));
          }

        /* Compute the line, rounded appropriately: */
        if (slope >= 0)
          { return yref + (slope*(x - xref)); }
        else
          { return yref - (slope*(xref - x)); }
      }
  }  

Float flt_interp_hi(Float xa, Float ya, Float xb, Float yb, Float x)
  { demand (xa < xb, "data points out of order"); 
    if (x == xa) 
      { return(ya); }
    else if (x == xb)
      { return(yb); }
    else 
      { ROUND_UP;
        /* Choose a reference point {(xref,yref)}: */
        Float xref, yref;
        if (x <= xa)
          { xref = xa; yref = ya; }
        else if (x >= xb)
          { xref = xb; yref = yb; }
        else if (FABS(x - xa) < FABS(x - xb))
          { xref = xa; yref = ya; }
        else 
          { xref = xb; yref = yb; }

        /* Compute the slope, rounded appropriately: */
        Float slope;
        if (x <= xref) 
          { /* We must round the slope down: */
            slope = -((ya - yb)/(xb - xa));
          }
        else
          { /* We must round the slope up: */
            slope = (yb - ya)/(-(xa - xb));
          }

        /* Compute the line, rounded appropriately: */
        if (slope >= 0)
          { return yref + (slope*(x - xref)); }
        else
          { return yref - (slope*(xref - x)); }
      }
  }  
  
Float flt_mix_lo(Float xa, Float wa, Float xb, Float wb)
  { demand((wa >= 0) && (wb >= 0), "negative weights");
    ROUND_UP;
    Float sw = wa + wb;
    if (sw == Zero) { return xa; }
    ROUND_DOWN;
    Float xr;
    if (xa < xb)
      { Float hdx = ((Half*xb) + (Half*(-xa)))*(wb/sw);
        xr = (xa + hdx) + hdx;
      }
    else
      { Float hdx = ((Half*xa) + (Half*(-xb)))*(wb/sw);
        xr = (xb + hdx) + hdx;
      }
    affirm ((xr >= xa) && (xr <= xb), "bug");
    return xr;
  }
    
Float flt_mix_hi(Float xa, Float wa, Float xb, Float wb)
  { demand((wa >= 0) && (wb >= 0), "negative weights");
    ROUND_DOWN;
    Float sw = wa + wb;
    if (sw == Zero) { return (xa); }
    ROUND_UP;
    Float xr;
    if (xa < xb)
      { Float hdx = ((Half*xb) + (Half*(-xa)))*(wb/sw);
        xr = (xa + hdx) + hdx;
      }
    else
      { Float hdx = ((Half*xa) + (Half*(-xb)))*(wb/sw);
        xr = (xb + hdx) + hdx;
      }
    affirm ((xr >= xa) && (xr <= xb), "bug");
    return xr;
  }

