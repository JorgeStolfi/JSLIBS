/* See {sample_conv.h}. */
/* Last edited on 2024-10-25 22:24:22 by stolfi */

#define _GNU_SOURCE
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

#define BT709_ENCODE_GAMMA 0.450
#define BT709_ENCODE_A     4.500
#define BT709_ENCODE_R     0.018053968511
#define BT709_ENCODE_B     0.099296826809
  /* Parameters for luminance encoding according to ITU-R Recommendation BT.709. */

float sample_conv_encode_BT709(float Y)
  { if(isnan(Y) || isinf(Y)) { return Y; }
    float a = fabsf(Y);
    if ((a == 0) || (a == 1)) { return Y; }
    if (a <= BT709_ENCODE_R)
      { a = (float)(a * BT709_ENCODE_A); }
    else
      { a = (float)((1 + BT709_ENCODE_B)*pow(a, BT709_ENCODE_GAMMA) - BT709_ENCODE_B); }
    return (Y < 0 ? -a : a);
  }

#define BT709_DECODE_GAMMA 2.222222222222
#define BT709_DECODE_A     0.222222222222
#define BT709_DECODE_R     0.081242858299
#define BT709_DECODE_B     0.099296826809
  /* Parameters for luminance decoding according to ITU-R Recommendation BT.709. */

float sample_conv_decode_BT709(float V)
  { if(isnan(V) || isinf(V)) { return V; }
    float a = fabsf(V);
    if ((a == 0) || (a == 1)) { return V; }
    if (a <= BT709_DECODE_R)
      { a = (float)(a * BT709_DECODE_A); }
    else
      { a = (float)pow((a + BT709_DECODE_B)/(1 + BT709_DECODE_B), BT709_DECODE_GAMMA); }
    return (V < 0 ? -a : a);
  }

float sample_conv_gamma(float z, double gamma, double bias)
  { demand(bias >= 0, "negative bias not implemented yet");
    if(isnan(z) || isinf(z)) { return z; }
    if (gamma == 1) { return z; }
    float a = fabsf(z);
    if ((a == 0) || (a == 1)) { return z; }
    if (bias == 0)
      { a = (float)pow(a, gamma); }
    else
      { double sg = sqrt(gamma);
        double c = pow(bias, 1/sg);
        double d = pow(bias, sg);
        double u = a*(1 - c) + c;
        double v = pow(u, gamma);
        double w = (v - d)/(1 - d);
        a = (float)w;
        if (a < 0) { a = 0; }
      }
    return (z < 0 ? -a : a);
  }

float sample_conv_log(float u, double bias, double uref, double logBase)
  {
    if ((! isfinite(logBase)) || (logBase == 0)) { return NAN; }
    if ((! isfinite(bias)) || (bias < 0)) { return NAN; }
    if ((! isfinite(uref)) || (uref <= 0)) { return NAN; }
    if (isnan(u) || (u < 0)) { return NAN; }
    if (u == +INFINITY) { return +INFINITY; }
    if (bias != 0) { u = hypot(u, bias); }
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
    if (u == -INFINITY) { u = 0.0; } else { u = uref*exp(u*logBase); }
    if (u == +INFINITY) { return +INFINITY; }
    if (u < bias) { return NAN; }
    if (bias != 0) { u = sqrt(u*u - bias*bias); }
    assert(isfinite(u));
    return (float)u;
  }

float sample_conv_interp(float z, int np, double U[], double V[])
  { if(isnan(z) || isinf(z)) { return z; }
    if (np == 0) { return z; }
    float a = fabsf(z);
    /* Binary search consecutive indices {ilo,ihi} such that {U[ilo] <= a <= U[ihi]}. */
    /* Assume that {U[-1]==V[-1]==0} and {U[np]==V[np]==1}. */
    int ilo = -1, ihi = np;
    if (a <= 0) 
      { ihi = 0; }
    else if (a >= 1)
      { ilo = np-1; }
    else
      { while (ihi - ilo >= 2)
          { int imd = (ilo + ihi)/2; 
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
  ( int iChan,           /* Channel index in input image. */
    int oChan,           /* Channel index in output image. */
    sample_uint32_t imin,   /* Minimum integer sample seen. */
    sample_uint32_t imax,   /* Maximum integer sample seen. */
    sample_uint32_t maxval, /* Maximum possible integer sample. */
    double lo,           /* Low end of float scaling range. */
    double hi,           /* High end of float scaling range. */
    float vmin,          /* Minimum float sample seen. */
    float vmax           /* Maximum float sample seen. */
  )
  { fprintf(stderr, "  converted int channel %d to float channel %d:\n", iChan, oChan);
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
    int *clo,
    int *chi,
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
  ( int iChan,           /* Channel index in input image. */
    int oChan,           /* Channel index in output image. */
    float vmin,          /* Minimum float sample seen. */
    float vmax,          /* Maximum float sample seen. */
    double lo,           /* Low end of float scaling range. */
    double hi,           /* High end of float scaling range. */
    int clo,             /* Number of samples seen below {lo}. */
    int chi,             /* Number of samples seen above {hi}. */
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
    for (int k = 0; k < chns; k++)
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
        for (int k = 0; k < chns; k++)
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
