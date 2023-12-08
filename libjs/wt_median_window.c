/* See wt_median_window.h */
/* Last edited on 2023-11-23 16:02:10 by stolfi */

#define wt_median_window_C_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <jsmath.h>
#include <wt_median.h>

#include <wt_median_window.h>

double wt_median_window
  ( int32_t nx,
    double x[],
    int32_t ix,
    int32_t nw,
    int32_t w[],
    bool_t interp,
    int32_t nk,
    int32_t kx[],
    int32_t *ns_P,
    double xs[],
    int32_t ws[]
  )
  {
    bool_t debug = FALSE;
    
    int32_t np = wt_median_window_index_set_update(nx, ix, nw, nk, kx);
    wt_median_index_set_sort(nx, x, nw, kx, np);
    int32_t ns = wt_median_gather_samples(nx, x, ix, nw, w, kx, xs, ws);
    if (debug)
      { fprintf(stderr, "    xs[0..%d] = ", ns-1);
        for (int32_t js = 0; js < ns; js++)
          { fprintf(stderr, " %24.16e", xs[js]);
            if (js > 0) { assert(xs[js] > xs[js-1]); }
          }
        fprintf(stderr, "\n");
      }
            
    double xm = wt_median_sorted(ns, xs, ws, interp);
    
    (*ns_P) = ns;
    return xm;
  }

int32_t wt_median_window_index_set_update
  ( int32_t nx,     /* Count of samples. */
    int32_t ix,     /* Index of first sample in window. */
    int32_t nw,     /* Window width. */
    int32_t nk,     /* (IN) Count of indices in {kx}. */
    int32_t kx[]    /* (IN/OUT) Sorted indices are {kx[0..nk-1]}. */
  )
  { bool_t debug = FALSE;
  
    demand(nw >= 0, "invalid window width {nw}");
    int32_t jx = ix + nw - 1; /* Last sample index in window. */
    demand((ix >= 0) && (jx < nx), "window spills outside of {0..nx-1}");
    
    demand((nk >= 0) && (nk <= nx), "invalid index set size {nk}");
    int32_t kx_min = INT32_MAX, kx_max = INT32_MIN;  /* Min and max valid indices in {kx[0..nk-1]}. */
    if (nk == 0)
      { if (debug) { fprintf(stderr, "  index set was empty\n"); } }
    else
      { /* Remove from {kx[0..nk-1]} the indices that are outside the window: */
        /* Also find min and max indices {kx_min,kx_max} in remaining set: */
        int32_t nk_new = 0;
        for (int32_t j = 0; j < nk; j++)
          { int32_t kxj = kx[j];
            demand((kxj >= 0) && (kxj < nx), "invalid index in {kx[..]}");
            if ((kxj >= ix) && (kxj <= jx))
              { kx[nk_new] = kxj; nk_new++;
                if (kxj < kx_min) { kx_min = kxj; }
                if (kxj > kx_max) { kx_max = kxj; }
              }
          }
        if (debug) { fprintf(stderr, "  removed %d indices, left %d", nk - nk_new, nk_new); }
        if (nk_new != 0)
          { if (debug) { fprintf(stderr, "  = {%d..%d}", kx_min, kx_max); }
            demand(nk_new == kx_max - kx_min + 1, "indices {kx[..]} are not consecutive");
          }
        if (debug) { fprintf(stderr, "\n"); }
        demand(nk_new <= nw, "there are repeated indices in {kx[..]}");
        nk = nk_new;
      }
    assert(nk <= nw);
    int32_t np = nw - nk;
    if (nk == 0)
      { for (int32_t i = 0; i < nw; i++) { kx[i] = ix + i; } }
    else
      { /* Append window indices  less than {kx_min} or greater than {kx_max}: */
        assert(kx_min >= 0);
        if (debug && (kx_min > ix)) { fprintf(stderr, "  adding {%d..%d}\n", ix, kx_min-1); }
        while(kx_min > ix) { kx_min--; kx[nk] = kx_min; nk++; }
        if (debug && (kx_max < jx)) { fprintf(stderr, "  adding {%d..%d}\n", kx_max+1, jx); }
        while(kx_max < jx) { kx_max++; kx[nk] = kx_max; nk++; }
        assert(nk == nw);
      }
    return np;
  }
