/* See {neuromat_filter_bandpass_log_gauss.h}. */
/* Last edited on 2024-01-06 08:15:25 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_filter_bandpass_log_gauss.h>

  
  

void neuromat_filter_bandpass_log_gauss_compute_parms
  ( double flo0,  double flo1,
    double fhi1,  double fhi0,
    double *fm_P, int32_t *np_P, double *sigma_P, double *mag_P,
    double *fsup_P,
    bool_t verbose
  )
  {
    demand((!isnan(flo0)) && (!isnan(flo1)) && (!isnan(fhi1)) && (!isnan(fhi0)), "invalid NAN arguments");
    demand((0 < flo0) && (flo0 < flo1) && (flo1 < fhi1) && (fhi1 < +INF), "invaild params sign/order");
    double tiny = 1.0e-6;
    
    /* Mean frequency of first hump: */
    double fm = flo1; /* Mean frequency of first hump. */
    
    /* Compute tentative {sigma}: */
    double sigma_hi = fabs(log(fhi0) - log(fhi1))/sqrt(-2*log(tiny));
    double sigma_lo = fabs(log(flo1) - log(flo0))/sqrt(-2*log(tiny));
    double sigma = fmin(sigma_hi, sigma_lo);
    if (verbose) { fprintf(stderr, "  %s: sigma = %16.12f\n", __FUNCTION__, sigma); }

    /* Compute number of bumps to cover {[flo1_fhi1]}: */
    double wd_aim = fabs(log(fhi1) - log(flo1));  /* Width of passband. */
    int32_t np = (int32_t)ceil(fmax(0, wd_aim/(2*sigma)) - 1.0e-6) + 1; /* Number of humps. */
    if (verbose) { fprintf(stderr, "  np = %d", np); }
    assert(np >= 1);
    demand(np < 100, "filter LG: band too broad compared to shoulder width");
    
    if (np > 1)
      { /* Adjust {sigma} to fit integer humps in passband: */
        double sigma_new = 0.5*wd_aim/(np - 1);
        affirm(sigma_new <= (1 + 1.0e-6)*sigma, "adjustment of {sigma} increased its value");
        sigma = sigma_new;
      }
    
    /* Decide the scaling factor {mag} based on number of bumps {np}: */
    double mag; 
    if (np == 1)
      { mag = 1.0; }
    else if ((np % 2) == 0)
      { double sum = 0.0;
        for (int32_t k = np/2; k >= 1; k--)
          { double d = 2*k - 1; /* Distance in multiples of sigma. */
            sum += 2*exp(-0.5*d*d);
          }
        mag = 1.0/sum;
      }
    else
      { double sum = 1.0;
        for (int32_t k = np/2; k >= 1; k--)
          { double d = 2*k; /* Distance in multiples of sigma. */
            sum += 2*exp(-0.5*d*d);
          }
        mag = 1.0/sum;
      }
    if (verbose) { fprintf(stderr, "  mag = %16.12f", mag); }
    
    double fsup = exp(log(fhi1) + 9*sigma);
    if (verbose) { fprintf(stderr, "  fsup = %16.12f\n", fsup); }

    (*fm_P) = fm;
    (*np_P) = np;
    (*sigma_P) = sigma;
    (*mag_P) = mag;
    (*fsup_P) = fsup;
  }
  
double neuromat_filter_bandpass_log_gauss_eval(double f, double fm, int32_t np, double sigma, double mag)
  { demand(fm > 0, "invalid {fm}");
    demand(np >= 0, "invalid {np}");
    demand((sigma != 0), "invalid {sigma}");
    f = fabs(f);
    if (f == 0) { return 0.0; }
    if (f == 0) { return 0; }
    double gtot = 0;
    for (uint32_t kp = 0;  kp < np; kp++)
      { double zf = fabs((log(f) - log(fm))/sigma - 2.0*kp);
        if (zf < 9.0) { double g = exp(-0.5*zf*zf); gtot += g; }
      }
    return mag*gtot;
  } 

