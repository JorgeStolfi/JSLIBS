/* See wt_median.h */
/* Last edited on 2024-11-23 06:10:09 by stolfi */

#define wt_median_C_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

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
  ( uint32_t ns,     /* Count of samples. */
    double xs[],     /* The samples are {xs[0..ns-1]}. */
    uint64_t ws[],     /* The respective weights are {ws[0..ns-1]}. */
    bool_t interp    /* Should interpolate the median? */
  )
  {
    bool_t debug = FALSE;
    
    /* Check weights and sample order and compute total weight: */
    demand(ns >= 0, "invalid number of samples");
    if (ns == 0) { return NAN; }    
    
    double xprev = -INF;
    int64_t wsum = 0;
    for (uint32_t i = 0;  i < ns; i++) 
      { uint64_t wsi = ws[i];
        demand(wsi <= wt_median_WEIGHT_SUM_MAX, "invalid weight");
        wsum += (int64_t)wsi; /* Should not overflow. */
        demand(wsum <= wt_median_WEIGHT_SUM_MAX, "sum of weights is too large");
        double xsi = xs[i];
        demand(isfinite(xsi), "samples must be finite");
        demand(xsi > xprev, "samples must be strictly increasing");
        xprev = xsi;
      }
    if (debug) { fprintf(stderr, "    wsum = %22ld\n", wsum); }
    assert(wsum > 0);
      
    /* Find indices {ia,ib} bracketing the median, and the corresponding {F} values: */
    int32_t ia, ib;
    int64_t Fa, Fb;
    if (ns == 1)
      { ia = ib = 0; Fa = Fb = 0; }
    else
      { /* Find {ia} = max {i} with {F(xs[i]) <= 0} : */
        ia = -1;
        Fa = -(int64_t)wsum;
        uint64_t wa_prev = 0;
        while (ia+1 < ns)
          { if (debug) { fprintf(stderr, "      ia = %d\n", ia); }
            if (debug) { fprintf(stderr, "      wa_prev = %22lu Fa =      %+22ld\n", wa_prev, Fa); }
            uint64_t wa_next = ws[ia+1];
            int64_t Fa_next = Fa + (int64_t)(wa_prev + wa_next);
            if (debug) { fprintf(stderr, "      wa_next = %22lu Fa_next = %+22ld\n", wa_next, Fa_next); }
            if (Fa_next > 0) { break; }
            ia++; Fa = Fa_next;
            if (Fa == 0) { break; }
            wa_prev = wa_next;
          }
        if (debug) { fprintf(stderr, "      ia final = %d\n", ia); }
        if (debug) { fprintf(stderr, "      Fa =      %+22ld\n", Fa); }
        assert((Fa <= 0) && (Fa >= -wsum));
        assert((ia >= 0) && (ia < ns));
        /* Find {ib} = min{i} with {F(xs[i]) >= 0}: */
        ib = (int32_t)ns;
        Fb = +wsum;
        uint64_t wb_prev = 0;
        while (ib > 0)
          { uint64_t wb_next = ws[ib-1];
            int64_t Fb_next = Fb - (int64_t)(wb_prev + wb_next);
            if (debug) { fprintf(stderr, "      ib = %d wb_prev = %22lu wb_next = %22lu Fb_next = %+22ld\n", ib, wb_prev, wb_next, Fb_next); }
            if (Fb_next < 0) { break; }
            ib--; Fb = Fb_next;
            if (Fb == 0) { break; }
            wb_prev = wb_next;
          }
        assert((Fb >= 0) && (Fb <= +wsum));
        assert((ib >= 0) && (ib < ns));
      }
    if (debug) { fprintf(stderr, "    ia = %d xs[ia] = %20.14f F(xs[ia]) = %+22ld\n", ia, xs[ia], Fa); }
    if (debug) { fprintf(stderr, "    ib = %d xs[ib] = %20.14f F(xs[ib]) = %+22ld\n", ib, xs[ib], Fb); }
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
        else if (labs(Fa) < labs(Fb))
          { return xa; }
        else if (labs(Fa) > labs(Fb))
          { return xb; }
        else
          { /* Choose the one with even rank: */
            return ((ia & 1) == 0 ? xa : xb);
          }
      }
  }

double wt_median_unsorted
  ( uint32_t n,       /* Count of samples. */
    double x[],       /* The samples are {x[0..n-1]}. */
    uint64_t w[],     /* The respective weights are {w[0..n-1]}. */
    bool_t interp,    /* Should interpolate the median? */
    uint32_t *ns_P,   /* (OUT) Number of distinct sample values with nonzero weight. */
    double xs[],      /* (WORK) Condensed sample table with {n} entries, or {NULL}. */
    uint64_t ws[],    /* (WORK) Condensed weight table with {n} entries, or {NULL}. */
    uint32_t kx[]     /* (WORK) Index table with {n} entries, or {NULL}. */
  )
  {
    bool_t debug = FALSE;
    demand (n >= 0, "invalid sample count {n}");

    uint32_t ns = 0;
    double xm = NAN;
    if (n > 0)
      { /* Create a trivial index set {kx[0..n-1]}: */
        for (uint32_t j = 0;  j < n; j++) { kx[j] = (uint32_t)j; }
        
        /* Sort the index set: */
        wt_median_index_set_sort(n, x, n, kx, 0);

        /* Gather distinct samples, sorted, with nonzero weights: */
        ns = wt_median_gather_samples(n, x, 0, n, w, kx, xs, ws);
        if (debug) { fprintf(stderr, "    ns = %d\n", ns); }
        if (debug)
          { fprintf(stderr, "    xs = \n"); 
            for (uint32_t ks = 0;  ks < ns; ks++) { fprintf(stderr, " %+14.8f", xs[ks]); }
            fprintf(stderr, "\n"); 
            fprintf(stderr, "    ws = \n"); 
            for (uint32_t ks = 0;  ks < ns; ks++) { fprintf(stderr, " %22lu", ws[ks]); }
            fprintf(stderr, "\n"); 
          }
    
        /* Compute median: */
        xm = wt_median_sorted(ns, xs, ws, interp);
      }
    /* Return results: */
    if (ns_P != NULL) { (*ns_P) = ns; }
    return xm;
  }

uint32_t wt_median_gather_samples
  ( uint32_t nx,     /* Count of samples. */
    double x[],      /* The samples are {x[0..nx-1]}. */
    uint32_t ix,     /* Lowest sample index in window. */
    uint32_t nw,     /* Window width. */
    uint64_t w[],    /* The window weights are {w[0..nw-1]}. */
    uint32_t kx[],   /* The previously sorted window indices are {kx[0..nw-1]}. */
    double xs[],     /* (OUT) the rearranged and condensed samples are {xs[0..*ns-1]}. */
    uint64_t ws[]    /* (OUT) the corresponding weights are {ws[0..*ns-1]}. */
  )
  {
    bool_t debug = FALSE;
    demand(ix + nw <= nx, "window {ix..ix+nw-1} overflows {0..nx-1}");
    demand(nw <= nx, "window too wide");
    
    double x_prev = -INF;
    uint32_t ns = 0;
    for (uint32_t ik = 0; ik < nw; ik++)
      { uint32_t jx = kx[ik];
        demand(jx < nx, "invalid index in {kx}");
        double xj = x[jx];
        demand(isfinite(xj), "sample value is not finite");
        demand(xj >= x_prev, "list {kx} is not sorted by sample value");
        demand((jx >= ix) && (jx < ix+nw), "index in {kx} lies outside the window");
        uint32_t iw = (uint32_t)(jx - ix);
        uint64_t wi = w[iw];
        demand(wi <= wt_median_WEIGHT_SUM_MAX, "invalid weight");
        if (wi > 0)
          { /* Append or condense to {xs,ws}: */
            if (xj != x_prev)
              { xs[ns] = xj; ws[ns] = 0; 
                ns++;
              }
            assert(ns > 0);
            ws[ns-1] += wi; /* Should not overflow.*/
            demand(ws[ns-1] <= wt_median_WEIGHT_SUM_MAX, "condensed weight is too large");
            if (debug) { fprintf(stderr, "    xj = %24.15e w[%d] = %22lu ns = %d ws[%d] = %22lu\n", xj, iw, wi, ns, ns-1, w[ns-1]); }
          }
        x_prev = xj;
      }
    return ns;
  }

void wt_median_index_set_sort(uint32_t nx, double x[], uint32_t nw, uint32_t kx[], uint32_t np)
  {
    demand(np <= nw, "invalid {np}");
    if (nw - np < 5)
      { wt_median_index_set_insertion_sort(nx, x, nw, kx); }
    else
      { wt_median_index_set_quick_sort(nx, x, nw, kx); }
  }
  
void wt_median_index_set_insertion_sort(uint32_t nx, double x[], uint32_t nw, uint32_t kx[])
  { for (int32_t i = 1;  i < nw; i++)
      { uint32_t txi = kx[i];
        int32_t j = i;
        while ((j > 0) && (x[kx[j-1]] > x[txi])) { kx[j] = kx[j-1]; j--; }
        kx[j] = txi;
      }
  }
  
void wt_median_index_set_quick_sort(uint32_t nx, double x[], uint32_t nw, uint32_t kx[])
  {
    auto int32_t compx(const void *a, const void *b);
    
    qsort(kx, nw, sizeof(uint32_t), compx);
    return;
    
    int32_t compx(const void *a, const void *b)
      { uint32_t *ia = (uint32_t *)a; assert((*ia) < nx);
        uint32_t *ib = (uint32_t *)b; assert((*ib) < nx);
        double xa = x[*ia], xb = x[*ib];
        if (xa < xb) 
          { return -1; }
        else if (xa > xb)
          { return +1; }
        else
          { return 0; }
      }
  }     

