/* See {sample_conv_hdyn.h}. */
/* Last edited on 2017-06-24 23:22:37 by stolfilocal */

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#include <sample_conv.h>
#include <sample_conv_hdyn.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#define BT_BIAS (sample_conv_BT709_BIAS) 
/* A {bias} parameter that approximates BT.709 when
  used with gamma {0.45} or {1/0.45}. */

interval_t sample_conv_hdyn_floatize
  ( sample_uint32_t iv,      /* Input integer sample ({Y} in formula). */ 
    double brght,          /* Brightness setting ({b}). */
    double ctrst,          /* Contrast setting ({c}). */ 
    double sigma,          /* Noise level ({s}). */
    double gamma,          /* Power law exponent ({g}). */ 
    /* !!! What about bias? !!! */
    sample_uint32_t black,   /* Black offset ({k}). */ 
    sample_uint32_t white,   /* White limit ({w}). */ 
    /* Statistical accumulators: */
    sample_uint32_t *imin,    /* Minimum {iv} seen. */
    sample_uint32_t *imax,    /* maximum {iv} seen. */
    int *clo,               /* Count of underexposed pixels. */
    int *chi,               /* Count of overexposed pixels. */
    float *vmin,            /* Minimum finite {LO(fv)} seen. */
    float *vmax             /* Maximum finite {HI(fv)} seen. */
  ) 
  { /* Update input statistics: */
    if ((imin != NULL) && (iv < (*imin))) { (*imin) = iv; }
    if ((imax != NULL) && (iv > (*imax))) { (*imax) = iv; }
    /* Floatize the input pixel value: */
    interval_t fv;   /* Result: */
    if (iv <= black)
      { LO(fv) = -INF; 
        HI(fv) = ((double)iv + 0.5 - black)/((double)white - black);
        if (clo != NULL) { (*clo)++; }
      }
    else if (iv >= white)
      { LO(fv) = ((double)iv - 0.5 - black)/((double)white - black); 
        HI(fv) = +INF;
        if (chi != NULL) { (*chi)++; }
      }
    else
      { LO(fv) = ((double)iv - 0.5 - black)/((double)white - black); 
        HI(fv) = ((double)iv + 0.5 - black)/((double)white - black);
      }
    double cwk = ctrst*(white - black);
    int d;
    for (d = 0; d < 2; d++)
      { if (isfinite(fv.end[d])) 
          { /* Apply gamma decoding: */
            if (gamma != 1.0) 
              { fv.end[d] = sample_conv_gamma((float)(fv.end[d]), gamma, BT_BIAS); }
            /* Undo brightness/contrast adjustment: */
            fv.end[d] = (fv.end[d] - brght)/cwk;
            /* Account for noise: */
            fv.end[d] += (d == 0 ? -3*sigma : +3*sigma);
          }
      }
    if ((vmin != NULL) && (isfinite(LO(fv))) && (LO(fv) < (*vmin))) { (*vmin) = (float)LO(fv); }
    if ((vmax != NULL) && (isfinite(HI(fv))) && (HI(fv) > (*vmax))) { (*vmax) = (float)HI(fv); }
    return fv;
  }

void sample_conv_hdyn_print_floatize_stats
  ( int iChan,           /* Channel index in input image. */
    int oChan,           /* Channel index in output image. */
    double brght,        /* Brightness setting ({b}). */
    double ctrst,        /* Contrast setting ({c}). */ 
    sample_uint32_t black, /* Black offset ({k}). */ 
    sample_uint32_t white, /* White limit ({w}). */ 
    sample_uint32_t imin,  /* Minimum integer sample seen. */
    sample_uint32_t imax,  /* Maximum integer sample seen. */
    int clo,             /* Count of underexposed pixels. */
    int chi,             /* Count of overexposed pixels. */
    float vmin,          /* Minimum float sample seen. */
    float vmax           /* Maximum float sample seen. */
  )
  { fprintf(stderr, "  converted int channel %d to float channel %d:\n", iChan, oChan);
    double fvmin = (0.0 - brght)/ctrst;
    double fvmax = (1.0 - brght)/ctrst;
    fprintf(stderr, "    mapped [ %u .. %u ] to [ %12.5e _ %12.5e ]\n", black, white, fvmin, fvmax);
    fprintf(stderr, "    actual input range  = [ %5u .. %5u ]\n", imin, imax);
    if ((clo > 0) || (chi > 0)) 
      { fprintf(stderr, "    underexposed pixels = %d  overexposed pixels = %d\n", clo, chi); }
    fprintf(stderr, "    actual finite output range = [ %12.5e _ %12.5e]\n", vmin, vmax);
  }

sample_uint32_t sample_conv_hdyn_quantize
  ( float fv, 
    double brght,          /* Brightness setting ({b}). */
    double ctrst,          /* Contrast setting ({c}). */ 
    double sigma,          /* Noise level ({s}). */
    double gamma,          /* Power law exponent ({g}). */ 
    sample_uint32_t black,   /* Black offset ({k}). */ 
    sample_uint32_t white,   /* White limit ({w}). */ 
    /* Statistical accumulators: */
    float *vmin,
    float *vmax, 
    int *clo,
    int *chi,
    sample_uint32_t *imin, 
    sample_uint32_t *imax
  )
  { double dfv; /* Value of {{fv} mapped for contrast, brightness, gamma. */
    if (isnan(fv)) 
      { dfv = 0.5; }
    else
      { if ((vmin != NULL) && (fv < (*vmin))) { (*vmin) = fv; }
        if ((vmax != NULL) && (fv > (*vmax))) { (*vmax) = fv; }
        /* Apply constrast and brightness: */
        dfv = fv * ctrst + brght;
        /* Clip to range: */
        if (dfv < 0.0) { dfv = 0.0;  if (clo != NULL) { (*clo)++; } }
        if (dfv > 1.0) { dfv = 1.0;  if (chi != NULL) { (*chi)++; } }
        /* Apply gamma encoding: */
        if (gamma != 1.0) 
          { dfv = sample_conv_gamma((float)dfv, 1/gamma, BT_BIAS);
          }
      }
    /* Quantize: */
    int64_t zv = (int64_t)black + (int64_t)floor(dfv*(double)((int64_t)white - (int64_t)black) + 0.5);
    demand((zv >= black) && (zv <= white), "bad {zv}");
    sample_uint32_t iv = (sample_uint32_t)zv;
    /* Update output range: */
    if ((imin != NULL) && (iv < (*imin))) { (*imin) = iv; }
    if ((imax != NULL) && (iv > (*imax))) { (*imax) = iv; }
    return iv;
  }
  
void sample_conv_hdyn_print_quantize_stats
  ( int iChan,           /* Channel index in input image. */
    int oChan,           /* Channel index in output image. */
    double brght,        /* Brightness setting ({b}). */
    double ctrst,        /* Contrast setting ({c}). */ 
    sample_uint32_t black, /* Black offset ({k}). */ 
    sample_uint32_t white, /* White limit ({w}). */ 
    float vmin,          /* Minimum float sample seen. */
    float vmax,          /* Maximum float sample seen. */
    int clo,             /* Number of samples seen below {lo}. */
    int chi,             /* Number of samples seen above {hi}. */
    sample_uint32_t imin,   /* Minimum integer sample seen. */
    sample_uint32_t imax    /* Maximum integer sample seen. */
  )
  { fprintf(stderr, "  converted float channel %d to int channel %d:\n", iChan, oChan);
    double fvmin = (0.0 - brght)/ctrst;
    double fvmax = (1.0 - brght)/ctrst;
    fprintf(stderr, "    mapped [ %12.5e _ %12.5e ] to [ %u .. %u ]\n", fvmin, fvmax, black, white);
    fprintf(stderr, "    actual input range  = [ %12.5e _ %12.5e ]\n", vmin, vmax);
    if ((clo > 0) || (chi > 0)) 
      { fprintf(stderr, "    clipped pixels: too low = %d  too high = %d\n", clo, chi); }
    fprintf(stderr, "    actual output range = [ %5u .. %5u ]\n", imin, imax);
  }

interval_t sample_conv_hdyn_merge_intervals(int n, interval_t fv[])
  {
    interval_t rv;  /* Result. */
    double sum_wa = 0; /* Sum of weight times average. */
    double sum_w = 0;  /* Sum of weights. */
    int k;
    for (k = 0; k < n; k++)
      { interval_t *fvk = &(fv[k]);
        demand(! isnan(LO(*fvk)), "interval LO cannot be NAN");
        demand(! isnan(HI(*fvk)), "interval HI cannot be NAN");
        demand(LO(*fvk) <= HI(*fvk), "input interval cannot be empty");
        if ((LO(*fvk) <= -INF) || (HI(*fvk) >= +INF))
          { /* Infinite interval -- ignore for now... */ }
        else
          { /* Convert to average and deviation: */
            double avgk = interval_mid(fvk);
            double devk = interval_rad(fvk)/3;
            double vark = devk*devk;
            double wk = 1/vark;
            sum_wa += wk*avgk;
            sum_w += wk;
          }
      }
    if (sum_w == 0) 
      { /* No data to decide, compute the intersection of infinite intervals: */
        LO(rv) = -INF; HI(rv) = +INF;
        for (k = 0; k < n; k++)
          { interval_t *fvk = &(fv[k]);
            if (LO(*fvk) <= HI(*fvk))
              { if ((LO(*fvk) <= -INF) && (HI(*fvk) < HI(rv))) { HI(rv) = HI(*fvk); }
                if ((HI(*fvk) >= +INF) && (LO(*fvk) > LO(rv))) { LO(rv) = LO(*fvk); }
              }
          }
      }
    else
      { double avg = sum_wa/sum_w;
        double var = 1/sum_w;
        double dev = sqrt(var);

        rv = interval_from_mid_rad(avg, 3*dev);
      }
    return rv;
  }
