#define PROG_NAME "test_wt_median"
#define PROG_DESC "test of {wt_median.h} and {wt_median_window.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-22 21:19:23 by stolfi */ 
/* Created on 2012-03-04 by J. Stolfi, UNICAMP */

#define test_wt_median_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <wt_table.h>
#include <wt_table_gaussian.h>
#include <wt_table_binomial.h>
#include <wt_table_triangular.h>
#include <wt_table_hann.h>
#include <wt_table_quantize.h>

#include <wt_median.h>
#include <wt_median_window.h>

int32_t main(int32_t argn, char **argv);

void twtm_test_median_sorted(uint32_t ns_max, bool_t verbose);
  /* Tests the procedure {wt_median_sorted} with up to {ns_max} elements. */

void twtm_test_median_unsorted(uint32_t n_max, bool_t verbose);
  /* Tests the procedure {wt_median_unsorted} with up to {n_max} samples and weights. */

void twtm_test_gather_samples(uint32_t nx_max, uint32_t nw_max, bool_t verbose);
  /* Tests the procedure {wt_median_gather_samples} with up to {nx_max}
    samples and {nw_max} weights and sample indices.
    Also tests {wt_median_index_set_sort}, {wt_median_index_set_quick_sort},
    {wt_median_index_set_insertion_sort}. */
  
void test_wt_median_window_index_set_update(uint32_t nx_max, uint32_t nw_max, bool_t verbose);
  /* Tests the procedure {wt_median_window_index_set_update} with up to {nx_max}
    samples and {nw_max} weights and sample indices. */

void twtm_check_gather_samples
  ( uint32_t nx,
    double x[],
    uint32_t ix,
    uint32_t nw,
    uint64_t w[],
    uint32_t ns,
    double xs[],
    uint64_t ws[],
    uint32_t kx[],
    bool_t verbose
  );
  /* Checks whether the vectors {xs[0..ns-1]} and {ws[0..ns-1]} are proper condensations of 
    the samples {x[kx[0..nw-1]]}.  Assumes that {w[[kx[i]-ix]} is the weight of {x[kx[i]]}
    for {i} in {0..nw-1}, and that {kx[0..nw-1]} is a permutation of {ix..ix+nw-1]}. */
    
void twtm_check_median
  ( uint32_t ns,
    bool_t sorted,
    double xs[],
    uint64_t ws[],
    bool_t interp,
    double xm,
    bool_t verbose
  );
  /* Checks whether {vm} is the correct output of {wt_median_sorted(ns,xs,ws,interp)}. */
  
void twtm_check_index_set_update
  ( uint32_t nx,
    uint32_t nw0,
    uint32_t kx0[],
    uint32_t ix1,
    uint32_t nw1,
    uint32_t kx1[],
    uint32_t nkept,
    bool_t verbose
  );
  /* Assumes that {wt_median_window_index_set_update} was applied to
    the index set {kx0[0..nw0-1]} with parameters {nx, ix1, nw1} resulted
    in {kx1[0..nw1-1]} and returning {nkept}. Checks if
    that is what was expected. */

void twtm_prxsws(uint32_t ns, char *xname, double xs[], char *wname, uint64_t ws[]);  
  /* Prints to stderr samples {xs[0..ns-1]} and weights
    {ws[0..ns-1]}, with names {xname} and {wname} respectively. */

uint32_t twtm_choose_num_samples(uint32_t ktry, uint32_t nx_max);
  /* Chooses a number of samples {nx} in {1..nx_max}. */
    
uint32_t twtm_choose_window_size(uint32_t ktry, uint32_t nw_max, uint32_t nx);
  /* Chooses a number of weights (the window size) {nw} in {1..min(nx,nw_max)}. */

void twtm_choose_samples(uint32_t nx, double xmax, double x[]);
  /* Fills {x[0..nx-1]} with random sample values in {[-xmax _ +xmax]}.  */

void twtm_choose_distinct_samples(uint32_t nx, double xmax, double x[]);
  /* Fills {x[0..nx-1]} with random but distnct
    sample values in {[-xmax _ +xmax]}. */

void twtm_choose_sorted_samples(uint32_t nx, double xmax, double x[]);
  /* Fills {x[0..nx-1]} with random but distnct and strictly
    increasing sample values in {[-xmax _ +xmax]}. */

void twtm_choose_weights(uint32_t nw, uint64_t wmin, uint64_t wmax, uint64_t w[]);
  /* Fills {w[0..nw-1]} with random weights in the range {[wmin _ wmax]}.
    Requires {wmin < wmax}.  Many weights, but not all, will be equal
    to either {wmin} or {wmax}. */

uint32_t twtm_choose_indices(uint32_t nx, uint32_t nk, uint32_t kx[]);
  /* Chooses an integer {ix} in {0..nx-nk} and fills {kx[0..nk-1]} with
    the set of {nk} consecutive indices {ix .. ix+nk-1}
    in random order.  Returns the lowest index {ix}. 
    Requires {1 <= nk <= nx}. */

void twtm_prkxw
  ( char *xname,
    uint32_t nx,
    double x[],
    char *wname,
    uint32_t nw,
    uint64_t w[],
    char *kname,
    uint32_t kx[]
  );
  /* Prints to stderr the indices {kx[0..nw-1]}, the samples {x[kx[0..nw-1]]}, and the weights
    {w[0..nw-1]} on one line, with names {xname} and {wname}
    respectively. */

void twtm_prxw(char *xname, uint32_t nx, double x[], uint32_t ix, char *wname, uint32_t nw, uint64_t w[]);
  /* Prints to stderr the samples {x[ix..ix+nw-1]} and weights
    {w[0..nw-1]} on one line, with names {xname} and {wname}
    respectively. */

void twtm_prkx(char *xname, uint32_t nk, uint32_t kx[], uint32_t nkept);
  /* Prints to stderr the perm {kx[0..nk-1]} on one line.
    Prints a '|' before {kx[nkept]}, if {nkept < nk}. */

int32_t main (int32_t argc, char **argv)
  {
    uint32_t nx_max = 1000; /* Max number of original samples. */
    uint32_t nw_max = 100;  /* Max window size. */
    uint32_t nprocs = 4;
    for (uint32_t iproc = 0; iproc < nprocs; iproc++) 
      { uint32_t ntries = 200;
        for (uint32_t ktry = 0; ktry < 200; ktry++)
          { bool_t verbose = ((ktry < 10) && (ktry <= ntries - 5));
            
            uint32_t nx = twtm_choose_num_samples(ktry, nx_max);
            assert(nx <= nx_max);
            
            uint32_t nw = twtm_choose_window_size(ktry, nx, nw_max);
            assert(nw <= nw_max);
            assert(nw <= nx);
            
            switch (iproc)
              { case 0: 
                  { twtm_test_median_sorted(nx, verbose); }
                  break;
                case 1: 
                  { twtm_test_gather_samples(nx, nw, verbose); }
                  break;
                case 2: 
                  { twtm_test_median_unsorted(nx, verbose); }
                  break;
                case 3: 
                  { test_wt_median_window_index_set_update(nx, nw, verbose); }
                  break;
                default: fprintf(stderr, "!! no case %d\n", iproc);
              }
          }
      }
    return 0;
  }
  
void twtm_test_median_sorted(uint32_t ns, bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "> --- %s ---\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, "  ns = %d ---\n", ns); }

    /* Generate {ns} distinct samples, sorted: */
    double *xs = talloc(ns, double);
    twtm_choose_sorted_samples(ns, 999.0, xs);

    /* Generate {ns} positive weights: */
    uint64_t *ws = talloc(ns, uint64_t);
    uint64_t wmin = 1;
    uint64_t wmax = (uint64_t)(2*ns);
    twtm_choose_weights(ns, wmin, wmax, ws);
    
    if (verbose && (ns <= 10)) { twtm_prxsws(ns, "xs", xs, "ws", ws); }

    /* Compute median with and without interpolation: */
    for (uint32_t ki = 0;  ki < 2; ki++)
      { bool_t interp = (ki != 0);
        if (verbose) { fprintf(stderr, "  ...... interp = %c ......\n", "FT"[interp]); }
        double xm = wt_median_sorted(ns, xs, ws, interp);
        if (verbose) { fprintf(stderr, "    xm =   %+14.8f = %24.15e\n", xm, xm); }
        twtm_check_median(ns, TRUE, xs, ws, interp, xm, verbose);
        if (verbose) { fprintf(stderr, "\n"); }
      }
    free(xs);
    free(ws);
    if (verbose) { fprintf(stderr, "< --- %s ---\n", __FUNCTION__); }
  }
         
void twtm_test_gather_samples(uint32_t nx, uint32_t nw, bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "> --- %s ---\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, "  nx = %d  nw = %d\n", nx, nw); }

    bool_t debug = TRUE;

    /* Choose {nx} random samples: */
    double *x = talloc(nx, double);
    twtm_choose_samples(nx, 999.0, x);

    /* Choose {nw} random non-negative weights: */
    uint64_t wmin = 0;
    uint64_t wmax = (uint64_t)(2*nx);
    uint64_t *w = talloc(nw, uint64_t);
    twtm_choose_weights(nw, wmin, wmax, w);
    
    /* Index list: */
    uint32_t *kx = talloc(nw, uint32_t);
    
    /* Gathered samples and weights: */
    uint32_t ns_max = nw;
    double *xs = talloc(ns_max, double);
    uint64_t *ws = talloc(ns_max, uint64_t);
    
    for (uint32_t which_sort = 0; which_sort < 3; which_sort++)
      { 
        if (debug) { fprintf(stderr, "  ...... which_sort = %d ......\n",  which_sort); }
        
        /* Choose {nw} random indices in {0..nx-1}: */
        uint32_t ix = twtm_choose_indices(nx, nw, kx);

        /* Choose sort procedure: */
        if (which_sort == 0)
          { /* Toss a ramdom {nkept}.  Should work even if incorrect. */
            uint32_t nkept = uint32_abrandom(0, nw);
            if (verbose) { fprintf(stderr, "  sort (auto, nkept = %d)\n", nkept); }
            wt_median_index_set_sort(nx, x, nw, kx, nkept);
          }
        else if (which_sort == 1)
          { if (verbose) { fprintf(stderr, "  sort (quick)\n"); }
            wt_median_index_set_quick_sort(nx, x, nw, kx);
          }
        else if (which_sort == 2)
          { if (verbose) { fprintf(stderr, "  sort (insertion)\n"); }
            wt_median_index_set_insertion_sort(nx, x, nw, kx);
          }
        else
          { assert(FALSE); }
        if (verbose && (nw <= 10)) { twtm_prkxw("x", nx, x, "w", nw, w, "k", kx); }

        /* Call the procedure: */
        uint32_t ns = wt_median_gather_samples(nx, x, ix, nw, w, kx, xs, ws);
        if (verbose && (ns <= 10)) { twtm_prxsws(ns, "xs", xs, "ws", ws); }

        /* Check the result: */
        twtm_check_gather_samples(nx, x, ix, nw, w, ns, xs, ws, kx, verbose);

      }
    free(x);
    free(w);
    free(kx);
    free(xs);
    free(ws);
    if (verbose) { fprintf(stderr, "< --- %s ---\n", __FUNCTION__); }
  }

void twtm_test_median_unsorted(uint32_t n, bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "> --- %s ---\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, "  n = %d\n", n); }
    
    /* Choose {n} random samples: */
    double *x = talloc(n, double);
    twtm_choose_distinct_samples(n, 999.0, x);
    
    /* Generate {n}  non-negative weights: */
    uint64_t wmin = 0;
    uint64_t wmax = (uint64_t)(2*n);
    uint64_t *w = talloc(n, uint64_t);
    twtm_choose_weights(n, wmin, wmax, w);

    if (verbose && (n <= 10)) { twtm_prxw("x", n, x, 0, "w", n, w); }
    
    /* Work vectors: */
    uint32_t ns_max = n;
    double *xs = talloc(ns_max, double);
    uint64_t *ws = talloc(ns_max, uint64_t);
    uint32_t *kx = talloc(ns_max, uint32_t);

    /* Compute median with and without interpolation: */
    for (uint32_t ki = 0;  ki < 2; ki++)
      { bool_t interp = (ki != 0);
        if (verbose) { fprintf(stderr, "  interp = %c\n", "FT"[interp]); }

        /* Compute median: */
        uint32_t ns = UINT32_MAX;
        double xm = NAN;
        if (verbose) { fprintf(stderr, "  ...... interp = %c ......\n", "FT"[interp]); }
        xm = wt_median_unsorted(n, x, w, interp, &ns, xs, ws, kx);
        if (verbose) { fprintf(stderr, "    xm = %+14.8f = %+24.15e ns = %d\n", xm, xm, ns); }
        twtm_check_median(n, FALSE, x, w, interp, xm, verbose);
        if (verbose && (ns <= 10)) { twtm_prxsws(ns, "xs", xs, "ws", ws); }
        twtm_check_median(ns, TRUE, xs, ws, interp, xm, verbose);
        if (verbose) { fprintf(stderr, "\n"); }
      }
    free(x);
    free(w);
    free(xs);
    free(ws);
    free(kx);
    if (verbose) { fprintf(stderr, "< --- %s ---\n", __FUNCTION__); }
  }

void test_wt_median_window_index_set_update(uint32_t nx, uint32_t nw, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "> --- %s ---\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, "  nx = %d nw = %d\n", nx, nw); }

    /* Current index set and saved copy thereof: */
    uint32_t nk_max = nx;
    uint32_t *kx = talloc(nk_max, uint32_t);
    uint32_t *kx0 = talloc(nk_max, uint32_t);
    
    /* Fill {kx} with a subrange {ix0..jx0} {0..nx-1}) in random order: */
    uint32_t nw0 = nw;
    uint32_t ix0 = twtm_choose_indices(nx, nw0, kx);
    uint32_t jx0 = ix0 + nw0 - 1;
    assert(jx0 < nx);
    if (verbose && (nw <= 15)) { twtm_prkx("initial", nw0, kx, UINT32_MAX); }

    uint32_t nupd = 10; /* Number of updates to try. */
    for (uint32_t kupd = 0; kupd < nupd; kupd++)
      { if (verbose) { fprintf(stderr, "  ...... update %d ......\n", kupd); }
        /* Paranoia: */
        demand(nw0 <= nx, "window size {nw0} too big for sample vector");

        /* Save a copy {kx0[1..nws-1]} of initial {kx[0..nw0-1]}: */
        for (uint32_t ik = 0;  ik < nw0; ik++) { kx0[ik] = kx[ik]; }
        
        /* Pick a new sub-window {ix1..ix1+nw1-1} that possibly overlaps with {ix0..ix0+nw0-1}: */
        uint32_t nw1 = (uint32_t)((2*nw0*kupd)/(nupd-1));
        if (nw1 > nx) { nw1 = nx; }
        uint32_t ix1 = uint32_abrandom(0, (uint32_t)(nx - nw1));

        /* Update the index set: */
        uint32_t nkept = wt_median_window_index_set_update(nx, nw0, kx, ix1, nw1);
        if (verbose && (nw1 <= 15)) { twtm_prkx("updated", nw1, kx, nkept); }

        /* Check whether it worked: */
        twtm_check_index_set_update(nx, nw0, kx0, ix1, nw1, kx, nkept, verbose);
        
        /* Set up for next update: */
        nw0 = nw1;
      } 
    if (verbose) { fprintf(stderr, "< --- %s ----\n", __FUNCTION__); }
  }

void twtm_check_index_set_update
  ( uint32_t nx,
    uint32_t nw0,
    uint32_t kx0[],
    uint32_t ix1,
    uint32_t nw1,
    uint32_t kx1[],
    uint32_t nkept,
    bool_t verbose
  )
  { if (verbose) { fprintf(stderr, "  > --- %s ---\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, "  nx = %u  nw0 = %u\n", nx, nw0); }
    if (verbose) { fprintf(stderr, "  ix1 = %u  nw1 = %u  nkept = %d\n",ix1, nw1, nkept); }
    
    /* Set {given[jx]} to true iff {jx} appears in {kx0[0..nw0-1]}: */
    bool_t given[nx];
    for (uint32_t jx = 0;  jx < nx; jx++) { given[jx] = FALSE; }
    for (uint32_t ik = 0;  ik < nw0; ik++) 
      { uint32_t jxi = kx0[ik];
        demand((jxi >= 0) && (jxi < nx), "invalid index in {kx0}");
        demand(! given[jxi], "repeated indices in {kx0}");
        given[jxi] = TRUE;
      }

    /* Check {kx1[0..kn1-1]}: */
    bool_t seen[nx]; /* To check that {kx1[0..nw1-1]} are all distinct. */
    for (uint32_t jx = 0;  jx < nx; jx++) { seen[jx] = FALSE; }
    for (uint32_t ik = 0;  ik < nw1; ik++)
      { uint32_t jxi = kx1[ik];
        demand((jxi >= ix1) && (jxi <= ix1 + nw1 - 1), "{kx1[ik]} outside window index range");
        demand(! seen[jxi], "repated index in {kx1}");
        seen[jxi] = TRUE;
        if (ik < nkept)
          { /* The first {nkept} elements must come from the input set: */
            demand(given[jxi], "inconsistent {nkept} - new listed as old");
          }
        else
          { /* The last {nw1-nkept} must be new: */
            demand(! given[jxi], "inconsistent {nkept} - old listed as new");
          }
      }
    if (verbose) { fprintf(stderr, "  < --- %s ----\n", __FUNCTION__); }
  }

uint32_t twtm_choose_num_samples(uint32_t ktry, uint32_t nx_max)
  { bool_t debug = FALSE;
    /* Choose {nx} and fill {x}: */
    uint32_t nx = 1 + (ktry % nx_max);
    assert(nx <= nx_max);
    if (ktry > 5) { nx += uint32_abrandom(0, (uint32_t)(nx_max - nx)); }
    if (debug) { fprintf(stderr, "    nx = %d\n", nx); }
    assert(nx <= nx_max);
    return nx;
  }

void twtm_choose_samples(uint32_t nx, double xmax, double x[])
  { 
    for (uint32_t i = 0;  i < nx; i++) 
      { x[i] = dabrandom(-xmax,+xmax);
        assert(fabs(x[i]) <= xmax);
      }
  }

void twtm_choose_sorted_samples(uint32_t nx, double xmax, double x[])
  { 
    double xr = 0.999*xmax*(nx == 1 ? 1.0 : 1.0/(nx+1)); /* Max perturbation to each sample. */
    for (uint32_t i = 0;  i < nx; i++) 
      { double xc = (2*((double)i+1)/(nx+1) - 1.0)*xmax;
        x[i] = xc + dabrandom(-xr, +xr);
        assert(fabs(x[i]) <= xmax);
        if (i > 0) { assert(x[i] > x[i-1]); }
      }
  }

void twtm_choose_distinct_samples(uint32_t nx, double xmax, double x[])
  { 
    /* Choose {nx} *sorted* distinct samples: */
    twtm_choose_sorted_samples(nx,xmax, x);
    if (nx >= 2)
      { /* Randomly permute them: */
        for (uint32_t i = 0; i < nx; i++)
          { uint32_t k = uint32_abrandom(i, (uint32_t)(nx-1));
            if (k != i) { double t = x[k]; x[k] = x[i]; x[i] = t; }
          }
      }
  }

uint32_t twtm_choose_window_size(uint32_t ktry, uint32_t nw_max, uint32_t nx) 
  {
    bool_t debug = FALSE;
    /* Choose {nw} and fill {w}: */
    uint32_t m = (nw_max < nx ? nw_max : nx);
    uint32_t nw = 1 + (ktry % m);
    assert(nw <= m);
    if (ktry > 5) { nw += uint32_abrandom(0, (uint32_t)(m - nw)); }
    if (debug) { fprintf(stderr, "  nw = %d\n", nw); }
    assert((nw <= nw_max) && (nw <= nx));
    return nw;
  }

void twtm_choose_weights(uint32_t nw, uint64_t wmin, uint64_t wmax, uint64_t w[]) 
  { demand(wmin < wmax, "bad or trivial range {wmin_wmax}");
    uint64_t dw = (wmax - wmin)/100;
    if ((wmax-wmin >= 2) && (dw == 0)) { dw = 1; }
    for (uint32_t kw = 0; kw < nw; kw++)
      { uint32_t p = (kw % 5);
        if (p == 1)
          { w[kw] = wmin; }
        else if (p == 3)
          { w[kw] = wmax; }
        else
          { w[kw] = uint64_abrandom(wmin+dw, wmax-dw); }
      }
  }
        
uint32_t twtm_choose_indices(uint32_t nx, uint32_t nk, uint32_t kx[])
  {
    bool_t debug = FALSE;
    
    demand((nk >= 1) && (nk <= nx), "invalid {nx,nk}");
    uint32_t jx0 = (nx - nk)/3;
    uint32_t jx1 = jx0 + nk - 1;
    if (debug) { fprintf(stderr, "  jx0 = %d jx1 = %d\n", jx0, jx1); }
    assert(jx0 <= jx1);
    if (nk == 1)
      { kx[0] = jx0; }
    else
      { for (uint32_t ik = 0; ik < nk; ik++)
          { /* Generate {kx[ik] By interpolation: */
            kx[ik] = jx0 + ik; 
            /* Shuffle it: */
            uint32_t rk = uint32_abrandom(0, ik);
            uint32_t t = kx[ik]; kx[ik] = kx[rk]; kx[rk] = t;
          }
      }
    return jx0;
  }

void twtm_check_median
  ( uint32_t ns,
    bool_t sorted,
    double xs[],
    uint64_t ws[],
    bool_t interp,
    double xm,
    bool_t verbose
  )
  { 
    bool_t debug = FALSE;
    if (verbose) { fprintf(stderr, "  > --- %s ---\n", __FUNCTION__); }
    if (verbose) { fprintf(stderr, "    ns = %d sorted = %c interp = %c\n", ns, "FT"[sorted], "FT"[interp]); }
    if (sorted)
      { /* Check strict order of {xs[0..ns-1]} and positive weights: */
        for (uint32_t ks = 0;  ks < ns; ks++)
          { demand(ws[ks] > 0, "invalid weights");
            if (ks > 0) { demand(xs[ks] > xs[ks-1], "samples not sorted"); }
          }
      }

    /* Compute weight sums and samples with best {F} on each side of zero: */
    bool_t is_a_sample = FALSE; /* True if {xm} is one of the {x[i]}. */
    int64_t Stot = 0; /* Sum of all weights. */
    int64_t Seq = 0;  /* Weight of {xm} if present. */
    
    /* Lower bound value {slo}, its condensed weight {wlo}. and its below-sum {Slo}. */
    double xlo = -INF; uint64_t wlo = 0; int64_t Slo = 0; 
    
    /* Upper bound value {shi}, its condensed weight {whi}. and its above-sum {Shi}. */
    double xhi = +INF; uint64_t whi = 0; int64_t Shi = 0;
    
    for (uint32_t ks = 0;  ks < ns; ks++)
      { int64_t Sall = Slo + Seq + Shi;
        if (debug) 
          { fprintf(stderr, "    ...... iteration %d ......\n", ks);
            fprintf(stderr, "    Stot = %+22ld\n", Stot); 
            fprintf(stderr, "    Seq =  %+22ld\n", Seq); 
            fprintf(stderr, "    is_a_sample = %c\n", "FT"[is_a_sample]); 
            fprintf(stderr, "    xlo = %+14.8f = %+24.15e  xhi = %+14.8f = %+24.15e\n", xlo, xlo, xhi, xhi); 
            fprintf(stderr, "    wlo = %22lu  whi = %22lu\n", wlo, whi); 
            fprintf(stderr, "    Slo = %+22ld  Shi = %+22ld\n", Slo, Shi);
            fprintf(stderr, "    Slo + Seq + Shi = %+22ld", Sall);
            fprintf(stderr, "\n");
          }
        assert(Sall == Stot);
              
        double xsk = xs[ks];
        uint64_t wsk = ws[ks];
        demand(wsk >= 0, "invalid (negative) weights");
        if (sorted) 
          { demand(wsk != 0, "invalid (zero) weights");
            if (ks > 0) { demand(xs[ks] > xs[ks-1], "samples not sorted"); }
          }
        if (wsk > 0)
          { Stot += (int64_t)wsk;
            if (xsk < xm)
              { Slo += (int64_t)wsk;
                if (xsk >= xlo) 
                  { if (xsk == xlo)
                      { if (sorted) { assert(FALSE); }
                        wlo += wsk;
                      } 
                    else
                      { xlo = xsk; wlo = wsk; }
                  }
              }
            else if (xsk > xm)
              { Shi += (int64_t)wsk;
                if (xsk <= xhi) 
                  { if (xsk == xhi)
                      { if (sorted) { assert(FALSE); }
                        whi += wsk;
                      } 
                    else
                      { xhi = xsk; whi = wsk; }
                  }
               }
             else
               { Seq += (int64_t)wsk; 
                 is_a_sample = TRUE;
               }
           }
      }

    assert((xlo < xm) && (xm < xhi));
    int64_t Fm = Slo - Shi;
    int64_t Flo = Fm - Seq - (int64_t)wlo; /* Estimated {F(xlo)}. */
    int64_t Fhi = Fm + Seq + (int64_t)whi; /* Estimated {F(xhi)}. */
    if (debug) 
      { fprintf(stderr, "    Flo = %+22ld\n", Flo);
        fprintf(stderr, "    Fm =  %+22ld\n", Fm);
        fprintf(stderr, "    Fhi = %+22ld\n", Fhi);
      }
    if (interp)
      { if (! is_a_sample)
          { demand ((Flo < 0) && (Fhi > 0), "wrong gap for median");
            /* Check interpolation: */
            assert(Seq == 0);
            assert(isfinite(xlo) && isfinite(xhi));
            assert((Flo < 0) && (Fhi > 0));
            double f = ((double)(0 - Flo))/((double)(Fhi - Flo));
            double xm_exp = (1-f)*xlo + f*xhi;
            if (verbose) { fprintf(stderr, "    f = %14.8f = %24.15e  xm_exp = %14.8f = %24.15e\n", f, f, xm_exp, xm_exp); }
            demand(fabs(xm - xm_exp) < 1.0e-12*(fabs(xlo) + fabs(xhi)), "incorrect interpolation");
          }
        else
          { demand(Slo == Shi, "should have interpolated");
            assert(Fm == 0);
          }
      }
    else
      { demand(is_a_sample, "median should be one of the samples"); 
        demand ((labs(Fhi) >= labs(Fm)) && (labs(Flo) >= labs(Fm)), "not optimal choice for median");
      }
    if (verbose) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
  }
        
void twtm_check_gather_samples
  ( uint32_t nx,
    double x[],
    uint32_t ix,
    uint32_t nw,
    uint64_t w[],
    uint32_t ns,
    double xs[],
    uint64_t ws[],
    uint32_t kx[],
    bool_t verbose
  )
  { 
    int64_t wsum_org = 0; /* Sum of all given weights. */
    for (uint32_t kw = 0;  kw < nw; kw++) { wsum_org += (int64_t)(w[kw]); }
    int64_t wsum_chk = 0; /* Sum of the condensed weights. */
    double xs_prev = -INF;
    for (uint32_t ks = 0;  ks < ns; ks++)
      { double xsk = xs[ks];
        demand(isfinite(xsk), "condensed sample is infinite or {NAN}");
        demand(xsk > xs_prev, "condensed samples out of order");
        uint64_t wsk_chk = 0;
        for (uint32_t ik = 0;  ik < nw; ik++)
          { uint32_t jx = kx[ik];
            demand((jx >= 0) && (jx < nx), "invalid index in {kx}");
            uint32_t iw = kx[ik] - ix;
            demand((iw >= 0) && (iw < nw), "index in {kx} out of window range");
            if (x[jx] == xsk) { wsk_chk += w[iw]; }
          }
        demand(wsk_chk == ws[ks], "condensed sample weight mismatch");
        wsum_chk += (int64_t)(ws[ks]);
        xs_prev = xsk;
      }
    demand(wsum_chk == wsum_org, "lost weight");
  }

void twtm_prxsws(uint32_t ns, char *xname, double xs[], char *wname, uint64_t ws[])      
  { fprintf(stderr, "\n");
    for (uint32_t k = 0;  k < ns; k++)
      { fprintf(stderr, "  %s[%4d] = %+22.14f = %+23.16e", xname, k, xs[k], xs[k]);
        fprintf(stderr, "  %s[%4d] = %22lu", wname, k, ws[k]);
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
  }
       
void twtm_prkxw(char *xname, uint32_t nx, double x[], char *wname, uint32_t nw, uint64_t w[], char *kname, uint32_t k[])      
  { fprintf(stderr, "  nx = %d  nw = %d\n", nx, nw);
    fprintf(stderr, "  %s[0..%d] =     ", kname, nw-1);
    for (uint32_t j = 0;  j < nw; j++) { fprintf(stderr, " %22u", k[j]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s[%s[0..%d]] =  ", xname, kname, nw-1);
    for (uint32_t j = 0;  j < nw; j++) { fprintf(stderr, " %+22.14f = %+24.16e", x[k[j]], x[k[j]]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s[0..%d] =     ", wname, nw-1);
    for (uint32_t j = 0;  j < nw; j++) { fprintf(stderr, " %22lu", w[j]); }
    fprintf(stderr, "\n");
  }
     
void twtm_prxw(char *xname, uint32_t nx, double x[], uint32_t ix, char *wname, uint32_t nw, uint64_t w[])      
  { uint32_t jx = ix+nw-1;
    fprintf(stderr, "  n = %d\n", nw);
    fprintf(stderr, "  %s[%d..%d] = ", xname, ix, jx);
    for (uint32_t k = 0; k < nw; k++) { fprintf(stderr, " %+22.14f", x[ix + k]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s[0..%d] = ", wname, nw-1);
    for (uint32_t k = 0;  k < nw; k++) { fprintf(stderr, " %22lu", w[k]); }
    fprintf(stderr, "\n");
  }
        
void twtm_prkx(char *xname, uint32_t nk, uint32_t kx[], uint32_t nkept)
  { fprintf(stderr, "  %s kx[0..%d] =", xname, nk-1); 
    for (uint32_t k = 0;  k < nk; k++) 
      { if (k == nkept) { fprintf(stderr, " |"); }
        fprintf(stderr, " %d", kx[k]);
      }
    if (nkept != UINT32_MAX) { fprintf(stderr, "  (nkept = %d)", nkept); }
    fprintf(stderr, "\n");
  }
