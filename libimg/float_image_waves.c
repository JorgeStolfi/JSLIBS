/* See {float_image_waves.h}. */

/* Last edited on 2023-01-17 15:50:58 by stolfi */ 
/* Created on 2009-06-02 by J. Stolfi, UNICAMP */

#define float_image_test_C_COPYRIGHT \
  "Copyright © 2009  by the State University of Campinas (UNICAMP)"

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
    double phase[],
    double squash
  )
  {
    /* Accumulate waves: */
    double res = 0.000;
    for (int32_t kf = 0; kf < NF; kf++)
      { /* Compute a longitudinal perturbation for the point: */
        double df = 1.0/(fx[kf]*fx[kf] + fy[kf]*fy[kf]);
        double dx = df*fx[kf];
        double dy = df*fy[kf];
        double rdang0 = 0.25*sin(0.5*M_PI*(fy[kf]*x - fx[kf]*y)); /* Transverse waviness. */
        double rdang1 = 0.40*rdang0*sin(0.5*M_PI*(fx[kf]*x + fy[kf]*y));
        double rdang = rdang0 + rdang1;
      
        /* Compute the wave's argument {ang}: */
        double ang = 2*M_PI*(fx[kf]*(x+rdang*dx) + fy[kf]*(y+rdang*dy) + phase[kf]); /* Phase of wave at {(x,y)}. */
        /* Compute {a}, the wave {kf}, and  {b = 1-a}, both in {(0 _ 1)}: */
        double a = amp[kf]*cos(ang); /* Wave amplitude, in {[0.001 _ 0.999]}. */
        res += a;
      }
      
    if (isfinite(squash) && (squash > 0))
      { /* Map {r} to {[-1 _ +1]}: */
        res = res/squash;
        res = res/hypot(1, res);
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
    double fmag = exp(log(0.50)/6); /* Factor for spatial freq and amp decrease. */
    
    /* Choose the amplitudes, phases, and frequency vectors: */
    double tt0 = dt/2; /* Angular direction of first wave. */
    double f0 = 0.40;
    double fxk = f0*cos(tt0); 
    double fyk = f0*sin(tt0);
    double ampk = 1.0;
    for (int32_t kf = 0; kf < NF; kf++)
      { /* Compute parameters of  and  {b = 1-a}, both in {(0 _ 1)}: */
        amp[kf] = ampk;
        fx[kf] = fxk;
        fy[kf] = fyk;
        phase[kf] = 0.15*(kf+1)*M_PI;
         /* Prepare for next wave: */
        double fxn = + fxk*cdt - fyk*sdt;
        double fyn = + fxk*sdt + fyk*cdt;
        fxk = fmag*fxn;
        fyk = fmag*fyn;
        ampk = fmag*ampk;
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
