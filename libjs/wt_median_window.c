/* See wt_median_window.h */
/* Last edited on 2024-11-22 21:19:08 by stolfi */

#define wt_median_window_C_COPYRIGHT \
  "Copyright � 2023  by the State University of Campinas (UNICAMP)"

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
  ( uint32_t nx,
    double x[],
    uint32_t ix,
    uint32_t nw,
    uint64_t w[],
    bool_t interp,
    uint32_t nk,
    uint32_t kx[],
    uint32_t *ns_P,
    double xs[],
    uint64_t ws[]
  )
  {
    bool_t debug = FALSE;
    demand(ix + nw <= nx, "window {ix..ix+nw-1} not contained in {0..nx-1}");
    demand(nk <= nx, "too many indices in current set");
    uint32_t nkept = wt_median_window_index_set_update(nx, nk, kx, ix, nw);
    wt_median_index_set_sort(nx, x, nw, kx, nkept);
    uint32_t ns = wt_median_gather_samples(nx, x, ix, nw, w, kx, xs, ws);
    if (debug)
      { fprintf(stderr, "    xs[0..%d] = ", ns-1);
        for (uint32_t js = 0;  js < ns; js++)
          { fprintf(stderr, " %24.16e", xs[js]);
            if (js > 0) { assert(xs[js] > xs[js-1]); }
          }
        fprintf(stderr, "\n");
      }
            
    double xm = wt_median_sorted(ns, xs, ws, interp);
    
    (*ns_P) = ns;
    return xm;
  }

uint32_t wt_median_window_index_set_update
  ( uint32_t nx,     /* Total count of samples. */
    uint32_t nk,     /* (IN/OUT) Count of indices in {kx}. */
    uint32_t kx[],   /* (IN/OUT) Sorted indices are {kx[0..nk-1]}. */
    uint32_t ix,     /* Index of first sample in new window. */
    uint32_t nw      /* Count of samples in new window. */
  )
  { bool_t debug = FALSE;
    demand(ix + nw <= nx, "window {ix..ix+nw-1} overflows {0..nx-1}");
    demand(nk <= nx, "too many indices in current set");
    uint32_t nkept; /* Count of indices that were preserved. */
    if (nw == 0)
      { if (debug) { fprintf(stderr, "  new index set is empty\n"); }
        nkept = 0;
      }
    else
      { demand((nk >= 0) && (nk <= nx), "invalid input index set size {nk}");
        uint32_t jx = (uint32_t)(ix + nw - 1); /* Last index in new window. */
        assert(jx < nx); /* Tested on entry. */
        uint32_t kx_min = (uint32_t)INT32_MAX;  /* Min retained elem of {kx[0..nk-1]}. */
        uint32_t kx_max = 0;  /* Max retained elem of {kx[0..nk-1]}. */
        uint32_t nk_new = 0;  /* Number of output indices defined so far. */
        if (nk == 0)
          { if (debug) { fprintf(stderr, "  input index set was empty\n"); } }
        else
          { /* Remove from {kx[0..nk-1]} the indices that are outside the new window: */
            /* Also find min and max indices {kx_min,kx_max} in remaining set: */
            for (uint32_t j = 0; j < nk; j++)
              { uint32_t kxj = kx[j];
                demand(kxj < nx, "invalid index in {kx[..]}");
                if ((kxj >= ix) && (kxj <= jx))
                  { /* Keep this element: */
                    kx[nk_new] = (uint32_t)kxj; nk_new++;
                    if (kxj < kx_min) { kx_min = kxj; }
                    if (kxj > kx_max) { kx_max = kxj; }
                  }
              }
            if (debug) 
              { fprintf(stderr, "  removed %d indices, left %d", nk - nk_new, nk_new);
                if (nk_new != 0) { fprintf(stderr, "  = {%d..%d}", kx_min, kx_max); }
                fprintf(stderr, "\n");
              }
            if (nk_new != 0) 
              { demand(nk_new == kx_max - kx_min + 1, "input indices were not consecutive"); }
          }
        demand(nk_new <= nw, "there were repeated indices in the input");
        nkept = nk_new; /* Count of new indices preserved. */
        if (debug) { fprintf(stderr, "  nkept = %u indices kept\n", nkept); }
        
        /* Add the new indices: */
        if (nk_new == 0)
          { /* Just store {ix..jx} into {kx[0..nw-1]}: */
            for (uint32_t i = 0; i < nw; i++) { kx[i] = ix + i; }
            nk_new = nw;
          }
        else
          { /* The kept indices must be a subrange {kx_min..kx_max} of {ix..jx}: */
            if (debug && (kx_min > ix)) { fprintf(stderr, "  adding {%d..%d}\n", ix, kx_min-1); }
            while(kx_min > ix) { kx_min--; kx[nk_new] = (uint32_t)kx_min; nk_new++; }
            if (debug && (kx_max < jx)) { fprintf(stderr, "  adding {%d..%d}\n", kx_max+1, jx); }
            while(kx_max < jx) { kx_max++; kx[nk_new] = (uint32_t)kx_max; nk_new++; }
          }
        assert(nk_new == nw);
      }
    return nkept;
  }
