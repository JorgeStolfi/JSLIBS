/* See {float_image_waves.h}. */

/* Last edited on 2023-01-10 20:07:36 by stolfi */ 
/* Created on 2009-06-02 by J. Stolfi, UNICAMP */

#define float_image_test_C_COPYRIGHT \
  "Copyright © 2009  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <math.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>

#include <float_image.h>

#include <float_image_waves.h>

double float_image_waves_eval
  ( double x,
    double y, 
    int32_t NF,
    double amp[], 
    double fx[], 
    double fy[], 
    double phase[]
  )
  {
    /* Accumulate waves: */
    double res = 0.000;
    for (int32_t kf = 0; kf < NF; kf++)
      { /* Compute {a}, the wave {kf} and  {b = 1-a}, both in {(0 _ 1)}: */
        double ph = 2*M_PI*(fx[kf]*x + fy[kf]*y + phase[kf]); /* Phase of wave at {(x,y)}. */
        double a = amp[kf]*cos(ph); /* Wave amplitude, in {[0.001 _ 0.999]}. */
        res += a;
      }
    return res;
  }

void float_image_waves_pick
  ( int32_t NF,
    double amp[], 
    double fx[], 
    double fy[], 
    double phase[],
    bool_t verbose
  )
  {
    double dt = 0.618*2*M_PI; /* Angle between directions component waves */
    double cdt = cos(dt);
    double sdt = sin(dt);
    double fmag = exp(log(1/3.0)/3); /* Factor for spatial freq and amp decrease. */
    
    /* Choose the amplitudes, phases, and frequency vectors: */
    double tt0 = dt/2; /* Angular direction of first wave. */
    double fxk = 0.33*cos(tt0); 
    double fyk = 0.33*sin(tt0);
    double ak = 1.0;
    for (int32_t kf = 0; kf < NF; kf++)
      { /* Compute {a}, the wave {kf} and  {b = 1-a}, both in {(0 _ 1)}: */
        amp[kf] = ak;
        fx[kf] = fxk;
        fy[kf] = fyk;
        phase[kf] = 0.15*(kf+1)*M_PI;
         /* Prepare for next wave: */
        double fxn = + fxk*cdt - fyk*sdt;
        double fyn = + fxk*sdt + fyk*cdt;
        fxk = fmag*fxn;
        fyk = fmag*fyn;
        ak = fmag*ak;
      }
    if (verbose)
      { fprintf(stderr, "wave parameters:\n");
        for (int32_t kf = 0; kf < NF; kf++) 
          { fprintf(stderr, "  %3d amp = %+15.12f", kf, amp[kf]);
            fprintf(stderr, "  fx = %+15.12f fy = %+15.12f", fx[kf], fy[kf]);
            fprintf(stderr, "  phase = %+15.12f\n", phase[kf]);
          }
      }
  }
