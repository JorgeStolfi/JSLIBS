/* See {neuromat_filter_tabulate_hartley_gains.h}. */
/* Last edited on 2024-01-03 14:26:06 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>
#include <neuromat_filter.h>

#include <neuromat_filter_tabulate_hartley_gains.h>
   
complex neuromat_filter_tabulate_hartley_gains_folded
  ( int32_t kf,
    int32_t nf,
    double fsmp,
    neuromat_filter_gain_t *gain,
    double fsup,
    bool_t verbose
  );
  /* If {fsup} is zero or negative, computes {W = gain(f)} where {f = kf*fsmp/nf}.
    If {fsup} is positive, computes the sum {W} of {gain(f + i*fsmp)}
    for all {i} such that such that {|f + i*fsmp| < fsup}. */

void neuromat_filter_tabulate_hartley_gains
  ( int32_t nf,
    double fsmp,
    neuromat_filter_gain_t *gain,
    double fsup,
    bool_t normalize,
    double H[],
    bool_t verbose
  )
  { 
    demand(fsmp > 0, "invalid {fsmp}");
    
    double Wmax = -INF;
    for (int32_t kf0 = 0; kf0 <= nf/2; kf0++)
      { complex w0 = neuromat_filter_tabulate_hartley_gains_folded(kf0, nf, fsmp, gain, fsup, verbose);
        Wmax = fmax(Wmax, cabs(w0));
        int32_t kf1 = (nf - kf0) % nf;  /* Index of the Hartley coeff with the same frequency. */
        if (kf0 == kf1)
          { /* Only one set of gains to add: */
            demand(fabs(cimag(w0)) < 1.0e-8, "Fourier coeff should be pure real");
            H[kf0] = creal(w0);
            if (verbose) 
              { fprintf(stderr, "  F[%d] = ( %16.12f + %16.12f * I )", kf0, creal(w0), cimag(w0)); 
                fprintf(stderr, " --> H[%d] = %16.12f\n\n", kf0, H[kf0]); 
              }
          }
        else
          { /* Two sets of gains to add: */
            complex w1 = neuromat_filter_tabulate_hartley_gains_folded(kf1, nf, fsmp, gain, fsup, verbose);
            Wmax = fmax(Wmax, cabs(w1));
            /* Force conjugation: */
            double wr = 0.5*(creal(w0) + creal(w1));
            double wi = 0.5*(cimag(w0) - cimag(w1));
            /* Convert Fourier coeff pair to Hartley coeff pair: */
            H[kf0] = wr + wi;
            H[kf1] = wr - wi;
            if (verbose) 
              { fprintf(stderr, "  F[%d] = ( %16.12f + %16.12f * I )", kf0, creal(w0), cimag(w0)); 
                fprintf(stderr, ", F[%d] = ( %16.12f + %16.12f * I )", kf1, creal(w1), cimag(w1)); 
                fprintf(stderr, " --> H[%d] = %16.12f, H[%d] = %16.12f\n\n", kf0, H[kf0], kf1, H[kf1]); 
              }
          }
      }
    if (normalize && (Wmax > 0) && (Wmax != 1.0))
      { /* Normalize to unit maximum gain: */
        if (verbose) { fprintf(stderr, "rescaling by 1 / %16.12f\n", Wmax); }
        for (int32_t kf = 0; kf < nf; kf++) { H[kf] /= Wmax; }
        Wmax = 1.0;
      }
    if (verbose) { fprintf(stderr, "max absolute gain is %16.12f\n", Wmax); }
  }
        
complex neuromat_filter_tabulate_hartley_gains_folded
  ( int32_t kf,
    int32_t nf,
    double fsmp,
    neuromat_filter_gain_t *gain,
    double fsup,
    bool_t verbose
  )
  {
    double f = (kf*fsmp)/nf;
    
    /* Compute range {imin..imax} for {i} such that {gain(f + i*fsmp)} should be added: */
    int32_t imin = (fsup <= 0 ? 0 : -(int32_t)floor((fsup + f)/fsmp - 1.0e-15));
    int32_t imax = (fsup <= 0 ? 0 : +(int32_t)floor((fsup - f)/fsmp + 1.0e-15));
    if (verbose && (imin != imax)) 
      { fprintf(stderr, "%6d  adding gain(f + i*fsmp) for f = %14.8f", kf, f);
        fprintf(stderr, "  i in {%d..%d}\n", imin, imax); 
      }

    /* Accumulate the complex gain {W = w(kf)}: */
    complex W = 0;
    for (int32_t i = imin; i <= imax; i++)
      { int32_t kfi = kf + i*nf;
        double fi = f + i*fsmp;
        if ((i == 0) || (fabs(fi) <= fsup)) 
          { complex gi = gain(fi);
            if (verbose && (imin != imax)) 
              { fprintf(stderr, "        i = %+3d  kf = %+6d  f = %+16.12f", i, kfi, fi);
                fprintf(stderr, "  g = ( %+16.12f + %+16.12f * I)\n", creal(gi), cimag(gi));
              }
            W += gi;
          }
      }
    if (verbose) { fprintf(stderr, "%6d  f = %16.12f  W = ( %+16.12f + %+16.12f * I)\n", kf, f, creal(W), cimag(W)); }
    if (verbose && (imin != imax)) { fprintf(stderr, "\n"); }
    return W;
  }
