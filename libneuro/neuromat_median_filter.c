/* See {neuromat_median_filter.h}. */
/* Last edited on 2024-11-18 12:36:26 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <rn.h>
#include <jsmath.h>
#include <wt_table.h>
#include <wt_table_generic.h>
#include <wt_median.h>
#include <wt_median_window.h>

#include <neuromat_median_filter.h>
 
void neuromat_median_filter_apply
  ( int32_t nt,
    int32_t ne, 
    double **val, 
    double **med, 
    int32_t nw,
    int32_t wt[],
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "  > %s nw = %d\n", __FUNCTION__, nw); }
    
    demand((nw >= 3) && ((nw % 2) == 1), "invalid window width {nw}");
    int32_t hw = nw/2; /* Window radius. */
   
    /* Raw samples and filtered samples: */
    int32_t nx = nt + 2*hw;    /* Number of samples to be filtered, including mirrored ends. */
    double *x = rn_alloc(nx);  /* Signal to be filtered, with mirrored ends. */
    double *s = rn_alloc(nt);  /* Median-smoothed signal. */
    
    /* Indices of window samples sorted by increasing value: */
    int32_t nk; 
    int32_t kx[nw]; /* Indices of previous window are {kx[0..nk-1]} in sorted {x} order. */
    
    /* Conslidated and sorted samples and weight: */
    double xs[nw];  /* Sorted distinct sample values. */
    int32_t ws[nw]; /* Consolidated weights of those values. */
    
    if (verbose) { fprintf(stderr, "    filtering electrodes \n"); }
    for (int32_t ie = 0; ie < ne; ie++)
      { if (verbose) { fprintf(stderr, "."); }
        /* Extract signal of electrode {ie} to {x[0..nt-1]}: */
        for (int32_t ix = 0; ix < nx; ix++)
          { int32_t it = ix - hw; /* Index of sample {x[jx]} in {val}. */
            /* Mirror about ends: */
            if (it < 0) { it = -it; }
            it = it % (2*nt - 2);
            if (it >= nt) { it = 2*nt-2 - it; }
            assert((it >= 0) && (it < nt));
            x[ix] = val[it][ie];
            demand(isfinite(x[ix]), "sample is infinite or {NAN}"); 
          }
            
        /* Apply running median filter: */
        nk = 0;
        for (int32_t it = 0; it < nt; it++)
          { /* Now {x[it+hw]} is the central sample of the window. */
            /* Check that the window fits in the sample vector: */
            assert(it+nw-1 < nx);
            int32_t ns = -1;
            s[it] = wt_median_window(nx, x, it, nw, wt, TRUE, nk, kx, &ns, xs, ws);
            assert(isfinite(s[it]));
            assert ((ns >= 1) && (ns <= nt));
            nk = nw;
          }
        /* Store the median signal into {med}: */
        for (int32_t it = 0; it < nt; it++) { med[it][ie] = s[it]; }
      }
    if (verbose) { fprintf(stderr, "\n"); }
    free(x);
    free(s);
    if (verbose) { fprintf(stderr, "  < %s\n", __FUNCTION__); }
    return;
  }
