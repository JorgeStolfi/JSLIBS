/* See {sample_conv.h}. */
/* Last edited on 2025-01-23 14:37:34 by stolfi */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <sample_conv.h>

#define BT709_ENC_EXPO  (0.45)
#define BT709_ENC_A     (4.5)
#define BT709_ENC_R     (0.018)
#define BT709_ENC_C     (0.099)
  /* Parameters for luminance encoding according to ITU-R Recommendation BT.709. */

float sample_conv_BT709_encode(float X)
  { if(isnan(X) || isinf(X)) { return X; }
    double a = fabsf(X);
    if ((a == 0) || (a == 1)) { return X; }
    if (a <= BT709_ENC_R)
      { a = (a * BT709_ENC_A); }
    else
      { a = ((1 + BT709_ENC_C)*pow(a, BT709_ENC_EXPO) - BT709_ENC_C); }
    return (float)(X < 0 ? -a : a);
  }

#define BT709_DEC_EXPO  (1/0.45)
#define BT709_DEC_A     (1/4.5)
#define BT709_DEC_R     (0.081)
#define BT709_DEC_C     (0.099) 
  /* Parameters for luminance decoding according to ITU-R Recommendation BT.709. */
  /* BT709_DEC_C (0.0992965876298) */

float sample_conv_BT709_decode(float V)
  { if(isnan(V) || isinf(V)) { return V; }
    double a = fabsf(V);
    if ((a == 0) || (a == 1)) { return V; }
    if (a <= BT709_DEC_R)
      { a = (a * BT709_DEC_A); }
    else
      { a = pow((a + BT709_DEC_C)/(1 + BT709_DEC_C), BT709_DEC_EXPO); }
    return (float)(V < 0 ? -a : a);
  }
  
/* SRGB */

#define sRGB_ENC_EXPO (1/2.4)
#define sRGB_ENC_A     (12.92)
#define sRGB_ENC_R     (0.0031308)
#define sRGB_ENC_C     (0.055)
  /* Parameters for linear component encoding according to IEC sRGB 1999 Amend 1. */

float sample_conv_sRGB_encode(float X)
  { if(isnan(X) || isinf(X)) { return X; }
    double a = fabsf(X);
    if ((a == 0) || (a == 1)) { return X; }
    if (a <= sRGB_ENC_R)
      { a = (a * sRGB_ENC_A); }
    else
      { a = ((1 + sRGB_ENC_C)*pow(a, sRGB_ENC_EXPO) - sRGB_ENC_C); }
    return (float)(X < 0 ? -a : a);
  }

#define sRGB_DEC_EXPO (2.4)
#define sRGB_DEC_A     (1/12.92)
#define sRGB_DEC_R     (0.04045)
#define sRGB_DEC_C     (0.055)
  /* Parameters for linear component decoding according to IEC sRGB 1999 Amend 1. */

float sample_conv_sRGB_decode(float V)
  { if(isnan(V) || isinf(V)) { return V; }
    double a = fabsf(V);
    if ((a == 0) || (a == 1)) { return V; }
    if (a <= sRGB_DEC_R)
      { a = (a * sRGB_DEC_A); }
    else
      { a = pow((a + sRGB_DEC_C)/(1 + sRGB_DEC_C), sRGB_DEC_EXPO); }
    return (float)(V < 0 ? -a : a);
  }

float sample_conv_log(float u, double bias, double uref, double logBase)
  {
    if ((! isfinite(logBase)) || (logBase == 0)) { return NAN; }
    if ((! isfinite(bias)) || (bias < 0)) { return NAN; }
    if ((! isfinite(uref)) || (uref <= 0)) { return NAN; }
    if (isnan(u) || (u < 0)) { return NAN; }
    if (u == +INFINITY) { return +INFINITY; }
    if (bias != 0) { u = (float)hypot(u, bias); }
    if (u == +INFINITY) { return +INFINITY; }
    if (u == 0.0) { return -INFINITY; }
    double ulog = log(u/uref)/logBase;
    return (float)ulog;
  }

float sample_conv_undo_log(float u, double bias, double uref, double logBase)
  {
    if ((! isfinite(logBase)) || (logBase == 0)) { return NAN; }
    if ((! isfinite(bias)) || (bias < 0)) { return NAN; }
    if ((! isfinite(uref)) || (uref <= 0)) { return NAN; }
    if (u == +INFINITY) { return +INFINITY; }
    if (isnan(u)) { return NAN; }
    if (u == -INFINITY) { u = 0.0; } else { u = (float)(uref*exp(u*logBase)); }
    if (u == +INFINITY) { return +INFINITY; }
    if (u < bias) { return NAN; }
    if (bias != 0) { u = (float)sqrt(u*u - bias*bias); }
    assert(isfinite(u));
    return (float)u;
  }

float sample_conv_interp(float z, int32_t np, double U[], double V[])
  { if(isnan(z) || isinf(z)) { return z; }
    if (np == 0) { return z; }
    float a = fabsf(z);
    /* Binary search consecutive indices {ilo,ihi} such that {U[ilo] <= a <= U[ihi]}. */
    /* Assume that {U[-1]==V[-1]==0} and {U[np]==V[np]==1}. */
    int32_t ilo = -1, ihi = np;
    if (a <= 0) 
      { ihi = 0; }
    else if (a >= 1)
      { ilo = np-1; }
    else
      { while (ihi - ilo >= 2)
          { int32_t imd = (ilo + ihi)/2; 
            /* assert((ilo < imd) && (imd < ihi)); */
            /* assert((0 <= imd) && (imd < np)); */
            if (a < U[imd])
              { ihi = imd; }
            else 
              { ilo = imd; }
          }
      }
    /* Get points {(ulo,vlo)} and {uhi,vhi)} that define the function for {a}: */
    double ulo, uhi, vlo, vhi;
    if (ilo < 0)   { ulo = vlo = 0; } else { ulo = U[ilo]; vlo = V[ilo]; }
    if (ihi >= np) { uhi = vhi = 1; } else { uhi = U[ihi]; vhi = V[ihi]; }
    
    /* Apply linear interpolation: */
    double s = (a - ulo)/(uhi - ulo); 
    float v = (float)((1-s)*vlo + s*vhi);
    return (z < 0 ? -v : +v);
  }

float sample_conv_floatize
  ( sample_uint32_t iv, 
    sample_uint32_t maxval, 
    bool_t isMask,
    double lo,
    double hi, 
    sample_uint32_t *imin, 
    sample_uint32_t *imax, 
    float *vmin,
    float *vmax
  )
  /* Convert integer intensity to float, kep stats: */
  { if ((imin != NULL) && (iv < (*imin))) { (*imin) = iv; }
    if ((imax != NULL) && (iv > (*imax))) { (*imax) = iv; }
    double rv;
    if (isMask) 
      { rv = ((double)iv)/((double)maxval); }
    else
      { rv = ((double)iv + 0.5)/((double)maxval + 1.0); }
    /* Convert to range {[lo_hi]}, making sure that 1.0 becomes {hi}: */
    float fv = (rv == 1.0 ? (float)hi : (float)(lo + rv*(hi - lo)));
    if ((vmin != NULL) && (fv < (*vmin))) { (*vmin) = fv; }
    if ((vmax != NULL) && (fv > (*vmax))) { (*vmax) = fv; }
    return fv;
  }

void sample_conv_print_floatize_stats
  ( int32_t iChan,           /* Channel index in input image. */
    int32_t oChan,           /* Channel index in output image. */
    sample_uint32_t imin,   /* Minimum integer sample seen. */
    sample_uint32_t imax,   /* Maximum integer sample seen. */
    sample_uint32_t maxval, /* Maximum possible integer sample. */
    double lo,           /* Low end of float scaling range. */
    double hi,           /* High end of float scaling range. */
    float vmin,          /* Minimum float sample seen. */
    float vmax           /* Maximum float sample seen. */
  )
  { fprintf(stderr, "  converted int32_t channel %d to float channel %d:\n", iChan, oChan);
    fprintf(stderr, "    mapped [ 0 .. %u ] to [ %12.5e _ %12.5e ]\n", maxval, lo, hi);
    fprintf(stderr, "    actual input range  = [ %5u .. %5u ]\n", imin, imax);
    fprintf(stderr, "    actual output range = [ %12.5e _ %12.5e]\n", vmin, vmax);
  }

sample_uint32_t sample_conv_quantize
  ( float fv, 
    sample_uint32_t maxval, 
    bool_t isMask,
    double lo,
    double hi, 
    float *vmin,
    float *vmax, 
    int32_t *clo,
    int32_t *chi,
    sample_uint32_t *imin, 
    sample_uint32_t *imax
  )
  { demand(! isnan(fv), "{fv} is NAN"); 
    demand(maxval > 0, "{maxval} is zero"); 
    if ((vmin != NULL) && (fv < (*vmin))) { (*vmin) = fv; }
    if ((vmax != NULL) && (fv > (*vmax))) { (*vmax) = fv; }
    sample_uint32_t smp;
    if (lo == hi)
      { smp = maxval/2; }
    else
      { double rv = (fv - lo)/(hi - lo);
        assert(! isnan(rv));
        if (rv <= 0.0) 
          { smp = 0;  if (clo != NULL) { (*clo)++; } }
        else if (rv >= 1)
          { smp = maxval; if (chi != NULL) { (*chi)++; } }
        else if (isMask) 
          { /* Map 0 to 0, 1 to {2*maxval}, linearly, with directed rounding: */ 
            double fsmp = rv*((double)maxval - 1.0e-12);
            fsmp = (rv < 0.5 ? ceil(fsmp - 0.5) : floor(fsmp + 0.5)); 
            assert((fsmp >= 0) && (fsmp <= maxval));
            smp = (sample_uint32_t)fsmp;
          }
        else
          { /* Round every interval {[iv/(maxval+1),(iv+1)/(maxval+1))} down to {iv}: */
            double fsmp = floor(rv*((double)maxval + 1.0 - 1.0e-12));
            assert((fsmp >= 0) && (fsmp <= maxval));
            smp = (sample_uint32_t)fsmp;
          }
      }
    if ((imin != NULL) && (smp < (*imin))) { (*imin) = smp; }
    if ((imax != NULL) && (smp > (*imax))) { (*imax) = smp; }
    return smp;
  }
  
void sample_conv_print_quantize_stats
  ( int32_t iChan,          /* Channel index in input image. */
    int32_t oChan,          /* Channel index in output image. */
    float vmin,             /* Minimum float sample seen. */
    float vmax,             /* Maximum float sample seen. */
    double lo,              /* Low end of float scaling range. */
    double hi,              /* High end of float scaling range. */
    int32_t clo,            /* Number of samples seen below {lo}. */
    int32_t chi,            /* Number of samples seen above {hi}. */
    sample_uint32_t maxval, /* Maximum possible integer sample. */
    sample_uint32_t imin,   /* Minimum integer sample seen. */
    sample_uint32_t imax    /* Maximum integer sample seen. */
  )
  { fprintf(stderr, "  converted float channel %d to int channel %d:\n", iChan, oChan);
    fprintf(stderr, "    mapped [ %12.5e _ %12.5e ] to [ 0 .. %u ]\n", lo, hi, maxval);
    fprintf(stderr, "    actual input range  = [ %12.5e _ %12.5e ]\n", vmin, vmax);
    fprintf(stderr, "    actual output range = [ %5u .. %5u ]\n", imin, imax);
    if ((clo > 0) || (chi > 0)) 
      { fprintf(stderr, "    clipped pixels: too low = %d  too high = %d\n", clo, chi); }
  }

void sample_conv_choose_maxval(uint32_t chns, sample_uint32_t imaxval[], sample_uint32_t maxmaxval, sample_uint32_t *omaxvalP)
  {
    sample_uint32_t imvmax = 0;
    sample_uint32_t imvmin = maxmaxval;
    for (int32_t k = 0; k < chns; k++)
      { sample_uint32_t imvk = imaxval[k];
        demand((imvk >= 1) && (imvk <= maxmaxval), "invalid channel maxval");
        if (imvmin > imaxval[k]) { imvmin = imaxval[k]; }
        if (imvmax < imaxval[k]) { imvmax = imaxval[k]; }
      }
    sample_uint32_t omaxval;
    if (imvmin == imvmax)
      { /* All channels have the same max value, use that: */
        omaxval = imvmax;
      }
    else if (imvmax == maxmaxval)
      { /* Some channel may have samples ranging up to {maxmaxval}, no other choice: */
        omaxval = maxmaxval;
      }
    else 
      { /* Pick the least common multiple, if it does not exceed {maxmaxval}: */
        omaxval = 1;
        for (int32_t k = 0; k < chns; k++)
          { sample_uint32_t imvk = imaxval[k];
            sample_uint32_t g = (sample_uint32_t)gcd(imvk, omaxval);
            sample_uint32_t d = imvk/g;
            /* Least common multiple so far would be {omaxval*d}: */
            if (maxmaxval/d < omaxval)
              { /* LCM would be too big: */
                omaxval = maxmaxval; 
                break;
              }
            else
              { /* Include {imvk} in the LCM: */
                omaxval = omaxval*d; 
                assert((omaxval % imvk) == 0);
              }
          }
       }
    (*omaxvalP) = omaxval;
  }
