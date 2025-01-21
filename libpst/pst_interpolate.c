/* See {pst_interpolate.h}  */

/* Created on 2011-06-17 by Jorge Stolfi, unicamp, <stolfi@dcc.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-01-16 12:56:41 by stolfi */
/* See the copyright and authorship notice at the end of this file.  */

#include <math.h>
#include <assert.h>
#include <stdint.h>

#include <affirm.h>
#include <float_image.h>
#include <pst_interpolate.h>

void pst_interpolate_two_values
  ( double v0, double w0,
    double v1, double w1,
    double *vRP, double *wRP
  )
  { *vRP = (v0+v1)/2;
    *wRP = 4.0/(1/w0 + 1/w1);
  }

void pst_interpolate_four_values
  (  double vm, double wm,
     double v0, double w0,
     double v1, double w1,
     double vp, double wp,
     double *vRP, double *wRP
  )
  {
    double dist_factor = 0.25;
    
    /* Extrapolate {vm,v0 --> va}: */
    double va = (3*v0 - vm)/2.0;
    double wa = dist_factor*4.0/(9.0/w0 + (1.0/(w0*wm)));
    
    /* Interpolate {v0,v1 --> vb}: */
    double vb = (v0+v1)/2.0;
    double wb = 4.0/(1.0/w0 + 1.0/w1);
    
    /* Extrapolate {v1,vp --> vc}: */
    double vc = (3*v1 - vp)/2.0;
    double wc = dist_factor*4.0/(9.0/w1 + (1.0/(w1*wp)));

    /* Weighted average of the estimates: */
    double vR, wR;
    if ((wa == 0) && (wb == 0) && (wc == 0))
      { if ((w0 != 0) || (w1 != 0)) 
          { /* Use the adjacent values: */
            vR = (w0*v0 + w1*v1)/(w0 + w1);
            wR = (w0 + w1)/2;
          }
        else
          { /* Give up: */ vR = wR = 0; }
      }
    else
      { /* Combine them: */
        wR = wa + wb + wc;
        vR = (wa*va + wb*vb + wc*vc)/(wa + wb + wc);
      }
    *wRP = wR; *vRP = vR; 
  }

