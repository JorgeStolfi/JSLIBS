/* See {neuromat_eeg_filter.h}. */
/* Last edited on 2023-12-14 09:01:40 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_eeg.h>
#include <neuromat_poly.h>
#include <neuromat_eeg_filter.h>

void neuromat_eeg_filter_apply
  ( int32_t nt, 
    int32_t ne, 
    double **val, 
    int32_t tdeg, 
    bool_t tkeep, 
    double G[],
    bool_t verbose
  )
  {
    /* Allocate the FFT work areas and plans: */
    double *in = (double*) fftw_malloc(sizeof(double)*(nt+4));
    double *out = (double*) fftw_malloc(sizeof(double)*(nt+4));
    fftw_plan pd = fftw_plan_r2r_1d(nt, in, out, FFTW_DHT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    fftw_plan pi = fftw_plan_r2r_1d(nt, out, in, FFTW_DHT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    
    /* Fitted polynomial coeffs and values: */
    double *P = (tdeg < 0 ? NULL : notnull(malloc(sizeof(double)*(tdeg+1)), "no mem"));
    double *s = (tdeg < 0 ? NULL : notnull(malloc(sizeof(double)*nt), "no mem"));
       
    for (int32_t ie = 0; ie < ne; ie++) 
      { /* Copy the signal {ie} into the FFT buffer, with mirrored apodized ends and zero padding: */
        for (int32_t it = 0; it < nt; it++) { in[it] = val[it][ie]; }
        
        if (tdeg >= 0)
          { /* Fit the polynomial {P}: */
            int32_t maxiter = 5;
            neuromat_poly_fit_robust(nt, NULL, in, NULL, maxiter, tdeg, P, NULL);
            if (verbose) 
              { fprintf(stderr, "  trend = ");
                for (int32_t r = 0; r <= tdeg; r++) { fprintf(stderr, " %+12.6f", P[r]); }
                fprintf(stderr, "\n");
              }
            /* Evaluate it: */
            neuromat_poly_eval_multi(tdeg, P, nt, NULL, s);
            /* Subtract it from the data: */
            for (int32_t it = 0; it < nt; it++) { in[it] -= s[it]; }
          }

        /* Transform to frequency domain: */
        fftw_execute(pd);

        /* Apply frequency filter: */
        for (int32_t kf0 = 0; kf0 <= nt-kf0; kf0++) 
          { int32_t kf1 = (nt - kf0) % nt; /* The other Hartley element with same absolute freq. */
            double h0 = out[kf0];
            double G0 = G[kf0];
            if (kf1 == kf0)
              { /* Coef {kf0} is the only one with that frequency: */
                out[kf0] = h0*G0;
              }
            else
              { /* Coefs {kf0,kf1} have the same frequency and get mixed: */
                double h1 = out[kf1]; 
                double G1 = G[kf1];
                /* Apply filter: */
                double Gp = G0 + G1;
                double Gm = G0 - G1;
                out[kf0] = 0.5*(h0*Gp + h1*Gm);
                out[kf1] = 0.5*(h1*Gp - h0*Gm);
              }
          }

        /* Return to time domain: */
        fftw_execute(pi);
        
        /* Return filtered signal & trend to {val}: */
        for (int32_t it = 0; it < nt; it++) 
          { /* Scale to preserve total power: */
            in[it] /= nt;
            /* Restore trend if any: */
            if ((tdeg >= 0) && tkeep) { in[it] += s[it]; }
            /* Save in {val}: */
            val[it][ie] = in[it];
          }
      }
      
    fftw_destroy_plan(pd);
    fftw_destroy_plan(pi);
    fftw_free(out);
    fftw_free(in);
    free(G);
    if (P != NULL) { free(P); }
    if (s != NULL) { free(s); }
    return;
  }
