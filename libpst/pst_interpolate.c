/* See {pst_interpolate.h}  */

/* Created on 2011-06-17 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-02-20 08:40:27 by stolfi */
/* See the copyright and authorship notice at the end of this file.  */

#include <math.h>
#include <assert.h>
#include <stdint.h>

#include <affirm.h>
#include <float_image.h>
#include <pst_interpolate.h>

void pst_interpolate_four_values
  (  double vm, double wm,
     double v0, double w0,
     double v1, double w1,
     double vp, double wp,
     double *vR_P, double *wR_P
  )
  {
    demand(isfinite(wm) && (wm >= 0), "invalid weight {wm}");
    demand(isfinite(w0) && (w0 >= 0), "invalid weight {w0}");
    demand(isfinite(w1) && (w1 >= 0), "invalid weight {w1}");
    demand(isfinite(wp) && (wp >= 0), "invalid weight {wp}");

    if (wm > 0) { demand(isfinite(vm), "invalid value {vm}"); }
    if (w0 > 0) { demand(isfinite(v0), "invalid value {v0}"); }
    if (w1 > 0) { demand(isfinite(v1), "invalid value {v1}"); }
    if (wp > 0) { demand(isfinite(vp), "invalid value {vp}"); }
    
    double vR, wR;
    if ((w0 == 0) || (w1 == 0))
      { vR = NAN; wR = 0.0; }
    else if ((wm == 0) && (wp == 0))
      { /* Linear interpolation {v0,v1}: */
        vR = (v0 + v1)/2;
        wR = 4/(1/w0 + 1/w1);
      }
    else if (wm == 0)
      { /* Quadratic interpolation {v0,v1,vp}: */
        vR = (3*v0 + 6*v1 - vp)/8; 
        wR = 64/(9/w0 + 36/w1 + 1/wp);
      }
    else if (wp == 0)
      { /* Quadratic interpolation {v1,v0,vm}: */
        vR = (3*v1 + 6*v0 - vm)/8; 
        wR = 64/(9/w1 + 36/w0 + 1/wm);
      }
    else
      { /* Cubic interpolation {vm,v0,v1,vp}: */
        vR = (-vm + 9*v0 + 9*v1 -vp)/16;
        wR = 256/(1/wm + 81/w0 + 81/w1 + 1/wp);
      }
      
    /* Check for overflow/underflow: */
    if ((! isfinite(vR)) || (! isfinite(wR)) || (wR == 0))
      { vR = NAN; wR = 0.0; }

    (*vR_P) = vR;
    (*wR_P) = wR; 
  }

