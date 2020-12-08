/* See {btc_price_series_local_avg.h} */
/* Last edited on 2015-04-20 22:15:12 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>

#include <bool.h>

#include <btc_price_series_local_fit_and_eval.h>

#include <btc_price_series_local_avg.h>

double btc_price_series_local_avg(int nd, double val[], int id, int hrad)
  { 
    bool_t debug = FALSE;
    
    if (hrad == 0) { return val[id]; }
    
    /* The procedure fits a polynomial by least squares to the
      logs of the values {val[kd]} in the window, weighted by a Hann-like window weight.
      Then it evaluates the fitted line
      at the current date.  Returns 0.0 if there is not enough
      data to fit a line. */
      
    /* Extract window values and compute windo weights: */
    int nw = 2*hrad + 1;
    double val_ex[nw]; /* Input value per day, within window. */
    double wht_ex[nw]; /* Hann weights, or 0 if no data. */
    int j;
    for (j = -hrad; j <= hrad; j++)
      { val_ex[hrad + j] = 0.0; /* By default. */
        wht_ex[hrad + j] = 0.0; /* By default. */
        int kd = id + j; 
        if ((kd >= 0) && (kd < nd))
          { double vj = val[kd];
            if (vj > 0)
              { double wj = 0.5*(1 + cos(M_PI*j/(hrad + 0.5)));
                val_ex[hrad + j] = vj;
                wht_ex[hrad + j] = wj;
                if (debug) { fprintf(stderr, "  %4d %18.5f %10.7f\n", j, val_ex[j], wht_ex[j]); }
              }
          }
      }
    
    /* Fit polynomial and evaluate at center: */
    double vsm = btc_price_series_local_fit_and_eval(1, hrad, val_ex, wht_ex);
    /* fprintf(stderr, " vsm = %18.5f\n", smo); */
    return vsm;
  }
          
