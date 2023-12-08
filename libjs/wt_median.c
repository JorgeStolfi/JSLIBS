/* See wt_median.h */
/* Last edited on 2023-11-23 15:53:40 by stolfi */

#define wt_median_C_COPYRIGHT \
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

double wt_median_sorted
  ( int32_t ns,     /* Count of samples. */
    double xs[],    /* The samples are {xs[0..ns-1]}. */
    int32_t ws[],   /* The respective weights are {ws[0..ns-1]}. */
    bool_t interp   /* Should interpolate the median? */
  )
  {
    bool_t debug = FALSE;
    
    /* Check weights and sample order and compute total weight: */
    demand(ns >= 0, "invalid number of samples");
    if (ns == 0) { return NAN; }    
    
    double xprev = -INF;
    int32_t wsum = 0;
    for (int32_t i = 0; i < ns; i++) 
      { int32_t wsi = ws[i];
        demand(wsi > 0, "invalid weight");
        demand(wsum <= wt_median_WSUM_MAX - wsi, "total weight too large");
        wsum += wsi;
        double xsi = xs[i];
        demand(isfinite(xsi), "samples must be finite");
        demand(xsi > xprev, "samples must be strictly increasing");
        xprev = xsi;
      }
    if (debug) { fprintf(stderr, "    wsum = %+d\n", wsum); }
    assert(wsum > 0);
      
    /* Find indices {ia,ib} bracketing the median, and the corresponding {F} values: */
    int32_t ia, ib;
    int32_t Fa, Fb;
    if (ns == 1)
      { ia = ib = 0; Fa = Fb = 0; }
    else
      { /* Find {ia} = max {i} with {F(xs[i]) <= 0} : */
        ia = -1;
        Fa = -wsum;
        int32_t wa_prev = 0;
        while (ia < ns-1)
          { int32_t wa_next = ws[ia+1];
            int32_t Fa_next = Fa + (wa_prev + wa_next);
            if (debug) { fprintf(stderr, "      ia = %d wa_prev = %d wa_prev = %d Fa_next = %+d\n", ia, wa_prev, wa_next, Fa_next); }
            if (Fa_next > 0) { break; }
            ia++; Fa = Fa_next;
            if (Fa == 0) { break; }
            wa_prev = wa_next;
          }
        assert((Fa <= 0) && (Fa >= -wsum));
        assert((ia >= 0) && (ia < ns));
        /* Find {ib} = min{i} with {F(xs[i]) >= 0}: */
        ib = ns;
        Fb = +wsum;
        int32_t wb_prev = 0;
        while (ib > 0)
          { int32_t wb_next = ws[ib-1];
            int32_t Fb_next = Fb - (wb_prev + wb_next);
            if (debug) { fprintf(stderr, "      ib = %d wb_prev = %d wb_prev = %d Fb_next = %+d\n", ib, wb_prev, wb_next, Fb_next); }
            if (Fb_next < 0) { break; }
            ib--; Fb = Fb_next;
            if (Fb == 0) { break; }
            wb_prev = wb_next;
          }
        assert((Fb >= 0) && (Fb <= +wsum));
        assert((ib >= 0) && (ib < ns));
      }
    if (debug) { fprintf(stderr, "    ia = %d xs[ia] = %20.14f F(xs[ia]) = %+d\n", ia, xs[ia], Fa); }
    if (debug) { fprintf(stderr, "    ib = %d xs[ib] = %20.14f F(xs[ib]) = %+d\n", ib, xs[ib], Fb); }
    assert((0 <= ia) && (ia <= ib) && (ib < ns));
    assert((Fa <= 0) && (Fb >= 0));
    if (ia == ib)
      { /* Median is exact and unique: */
        assert((Fa == 0) && (Fb == 0));
        return xs[ia];
      }
    else
      { /* Median is between two samples {xs[ia],xs[ib]}: */
        assert(ia+1 == ib);
        assert((Fa < 0) && (Fb > 0));
        double xa = xs[ia], xb = xs[ib];
        assert(xa < xb);
        if (interp)
          { /* Estimate {xm} where {F(xm)==0} by interpolation: */
            double f = ((double)Fb)/((double)(Fb - Fa));
            double xm = f*xs[ia] + (1-f)*xs[ib];
            if (xm <= xa) { xm = nextafter(xm, +INF); }
            if (xm >= xb) { xm = nextafter(xm, -INF); }
            return xm;
          }
        else if (abs(Fa) < abs(Fb))
          { return xa; }
        else if (abs(Fa) > abs(Fb))
          { return xb; }
        else
          { /* Choose the one with even rank: */
            return ((ia & 1) == 0 ? xa : xb);
          }
      }
  }

double wt_median_unsorted
  ( int32_t n,      /* Count of samples. */
    double x[],     /* The samples are {x[0..n-1]}. */
    int32_t w[],    /* The respective weights are {w[0..n-1]}. */
    bool_t interp,  /* Should interpolate the median? */
    int32_t *ns_P,  /* (OUT) Number of distinct sample values with nonzero weight. */
    double *xs,     /* (WORK) Condensed sample table with {n} entries, or {NULL}. */
    int32_t *ws,    /* (WORK) Condensed weight table with {n} entries, or {NULL}. */
    int32_t *kx     /* (WORK) Index table with {n} entries, or {NULL}. */
  )
  {
    bool_t debug = FALSE;
    
    demand (n >= 0, "invalid sample count {n}");
    int32_t ns = 0;
    double xm = NAN;
    if (n > 0)
      { /* Provide work arrays if needed, and remember those to be freed: */
        double *xs_alloc = NULL;
        if (xs == NULL) { xs_alloc = talloc(n, double); xs = xs_alloc; }
        int32_t *ws_alloc = NULL;
        if (ws == NULL) { ws_alloc = talloc(n, int32_t); ws = ws_alloc; }
        int32_t *kx_alloc = NULL;
        if (kx == NULL) { kx_alloc = talloc(n, int32_t); kx = kx_alloc; }

        /* Create a trivial index set {kx[0..n-1]}: */
        for (int32_t j = 0; j < n; j++) { kx[j] = j; }
        
        /* Sort the index set: */
        wt_median_index_set_sort(n, x, n, kx, n);

        /* Gather distinct samples, sorted, with nonzero weights: */
        ns = wt_median_gather_samples(n, x, 0, n, w, kx, xs, ws);
        if (debug) { fprintf(stderr, "    ns = %d\n", ns); }
        if (debug)
          { fprintf(stderr, "    xs = \n"); 
            for (int32_t ks = 0; ks < ns; ks++) { fprintf(stderr, " %+14.8f", xs[ks]); }
            fprintf(stderr, "\n"); 
            fprintf(stderr, "    ws = \n"); 
            for (int32_t ks = 0; ks < ns; ks++) { fprintf(stderr, " %14d", ws[ks]); }
            fprintf(stderr, "\n"); 
          }
    
        /* Compute median: */
        xm = wt_median_sorted(ns, xs, ws, interp);
        
        /* Free internally allocated storage: */
        if (xs_alloc != NULL) { free(xs_alloc); }
        if (ws_alloc != NULL) { free(ws_alloc); }
        if (kx_alloc != NULL) { free(kx_alloc); }
      }
    /* Return results: */
    if (ns_P != NULL) { (*ns_P) = ns; }
    return xm;
  }

int32_t wt_median_gather_samples
  ( int32_t nx,     /* Count of samples. */
    double x[],     /* The samples are {x[0..nx-1]}. */
    int32_t ix,     /* Lowest sample index in window. */
    int32_t nw,     /* Window width. */
    int32_t w[],    /* The window weights are {w[0..nw-1]}. */
    int32_t kx[],   /* The previously sorted window indices are {kx[0..nw-1]}. */
    double xs[],    /* (OUT) the rearrranged and condensed samples are {xs[0..*ns-1]}. */
    int32_t ws[]    /* (OUT) the corresponding weights are {ws[0..*ns-1]}. */
  )
  {
    bool_t debug = FALSE;
    
    double x_prev = -INF;
    int32_t ns = 0;
    for (int32_t ik = 0; ik < nw; ik++)
      { int32_t jx = kx[ik];
        demand((jx >= 0) && (jx < nx), "invalid index in {kx}");
        double xj = x[jx];
        demand(isfinite(xj), "sample value is not finite");
        demand(xj >= x_prev, "list {kx} is not sorted by sample value");
        int32_t iw = jx - ix;
        demand((iw >= 0) && (iw < nw), "index in {kx} lies outside the window");
        int32_t wi = w[iw];
        demand(wi >= 0, "weight is negative");
        if (wi > 0)
          { /* Append or condense to {xs,ws}: */
            if (xj != x_prev)
              { xs[ns] = xj; ws[ns] = 0; 
                ns++;
              }
            assert(ns > 0);
            demand(wi <= wt_median_WSUM_MAX - ws[ns-1], "condensed weight is too big");
            ws[ns-1] += wi;
            if (debug) { fprintf(stderr, "    xj = %24.15e wi = %d ns = %d ws[%d] = %d\n", xj, wi, ns, ns-1, w[ns-1]); }
          }
        x_prev = xj;
      }
    return ns;
  }

void wt_median_index_set_sort(int32_t nx, double x[], int32_t nw, int32_t kx[], int32_t np)
  {
    if (np < 5)
      { wt_median_index_set_insertion_sort(nx, x, nw, kx); }
    else
      { wt_median_index_set_quick_sort(nx, x, nw, kx); }
  }
  
void wt_median_index_set_insertion_sort(int32_t nx, double x[], int32_t nw, int32_t kx[])
  { for (int32_t i = 1; i < nw; i++)
      { int32_t txi = kx[i];
        int32_t j = i;
        while ((j > 0) && (x[kx[j-1]] > x[txi])) { kx[j] = kx[j-1]; j--; }
        kx[j] = txi;
      }
  }
  
void wt_median_index_set_quick_sort(int32_t nx, double x[], int32_t nw, int32_t kx[])
  {
    auto int32_t compx(const void *a, const void *b);
    
    qsort(kx, nw, sizeof(int32_t), compx);
    return;
    
    int32_t compx(const void *a, const void *b)
      { int32_t *ia = (int32_t *)a; assert((*ia) >= 0 && ((*ia) < nx));
        int32_t *ib = (int32_t *)b; assert((*ib) >= 0 && ((*ib) < nx));
        double xa = x[*ia], xb = x[*ib];
        if (xa < xb) 
          { return -1; }
        else if (xa > xb)
          { return +1; }
        else
          { return 0; }
      }
  }     

