/* See {neuromat_median_filter.h}. */
/* Last edited on 2023-11-04 02:38:28 by stolfi */

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
#include <wt_median.h>

#include <neuromat_median_filter.h>
 
void neuromat_median_filter_apply
  ( int32_t nt,
    int32_t ne, 
    double **val, 
    double **med, 
    int32_t nw,
    bool_t verbose
  )
  {
    demand((nw >= 3) && ((nw % 2) == 1), "invalid window width {nw}");
    int32_t hw = nw/2; /* Window radius. */
   
    /* The window weights are {wt[0..2*hw_cur]}, where {hw_cur} is in {0..hw}: */
    int32_t hw_cur = -1;  /* Current window radius. */ 
    int32_t nw_cur = -1;  /* Current window width. */ 
    double *wt = rn_alloc(nw);  /* Allocated with max size, use {wt[0..nw_cur-1]}. */
    
    double *x = rn_alloc(nt); /* Signal to be filtered. */
    double *s = rn_alloc(nt); /* Median-smoothed signal. */
    
    /* Indices of window samples sorted by increasing value: */
    int32_t nk; 
    int32_t kx[nw]; /* Indices of previous window are {kx[0..nk-1]} in sorted {x} order. */
    
    for (int32_t ie = 0; ie < ne; ie++)
      { /* Extract electrode {ie} signal to {x[0..nt-1]}: */
        for (int32_t it = 0; it < nt; it++) { x[it] = val[it][ie]; }
        
        /* Apply running median filter: */
        nk = 0;
        for (int32_t it = 0; it < nt; it++)
          { /* Compute the window radius {hwi} and width {nwi} to use for this sample: */
            int32_t hwi = (it < hw ? it : (it > nt - hw - 1 ?  nt - 1 - it : hw));
            int32_t nwi = 2*hwi + 1; 
            assert(nwi <= nw);
            /* Make sure that the weight table has size {2*hwi+1} */
            if (hw_cur != hwi)
              { /* Recompute table: */
                hw_cur = hwi; nw_cur = nwi;
                bool_t normalize = TRUE; /* Normalize table to unit weight. */
                wt_table_fill_hann(nw_cur, wt, normalize);
              }

            s[it] = wt_median(nt, x, it, nw_cur, wt, TRUE, &nk, kx);
          }
        /* Store the smoothed signal: */
        for (int32_t it = 0; it < nt; it++) { med[it][ie] = s[it]; }
      }
    free(x);
    free(s);
    return;
  }
