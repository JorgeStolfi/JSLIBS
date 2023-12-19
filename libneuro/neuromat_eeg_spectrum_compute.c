/* See {neuromat_eeg_spectrum_compute.h}. */
/* Last edited on 2023-12-14 08:35:11 by stolfi */

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
#include <neuromat_eeg_spectrum_compute.h>

double **neuromat_eeg_spectrum_compute(int32_t nt, int32_t ne, double **val, int32_t kfmax, bool_t verbose)  
  {
    demand(kfmax <= nt/2, "parameter {kfmax} exceeds the Nyquist frequency");
    
    /* Allocate the output: */
    double **pwr = notnull(malloc(ne*sizeof(double*)),"no mem"); 
    
    /* Allocate the work area: */
    double *in = (double*) fftw_malloc(sizeof(double) * nt);
    double *out = (double*) fftw_malloc(sizeof(double) * nt);
    
    /* Apodizing window: */
    double *W = notnull(malloc(nt*sizeof(double)), "no mem");
    for (int32_t it = 0; it < nt; it++) { W[it] = 0.5*(1 - cos(2*M_PI*(it+0.5)/nt)); }
    
    /* Do the row transforms: */
    fftw_plan px = fftw_plan_r2r_1d(nt, in, out, FFTW_DHT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    for (int32_t ie = 0; ie < ne; ie++)
      { /* Allocate and compute the spectrum of electrode {ie}: */
        pwr[ie] = notnull(malloc((kfmax+1)*sizeof(double)), "no mem"); 
        double sum2_in = 0;
        for (int32_t it = 0; it < nt; it++) 
          { double vt = val[it][ie] * W[it];
            sum2_in += vt*vt;
            in[it] =  vt;
          }
        double rms_in = sqrt(sum2_in / nt);
        fftw_execute(px);
        if (verbose) { fprintf(stderr, "  electrode %d\n", ie); }
        double sum2_out = 0;
        for (int32_t kf0 = 0; kf0 <= kfmax; kf0++) 
          { /* Get the index {kf1} of the coefficient with same actual frequency {kf0}: */
            int32_t kf1 = (nt - kf0) % nt;
            /* Compute the total power in that frequency: */
            double c0 = out[kf0];
            double c1 = out[kf1];
            if (verbose) 
              { fprintf(stderr, "    kf0 = %5d out[kf0] = %24.15e", kf0, c0); 
                if (kf1 != kf0) { fprintf(stderr, "  kf1 = %5d out[kf1] = %24.15e", kf1, c1); }
                fprintf(stderr, "\n");
              }
            double pf = (c0*c0 + (kf0 != kf1 ? c1*c1 : 0.0))/nt;
            sum2_out += pf;
            pwr[ie][kf0] = pf/nt;
          }
        double rms_out = sqrt(sum2_out/nt);
        if (verbose) { fprintf(stderr, "rms ratio out/in = %24.15e\n", rms_out/rms_in); }
      }
    fftw_destroy_plan(px);
    
    fftw_free(out);
    fftw_free(in);
    free(W);
    return pwr;
  }

