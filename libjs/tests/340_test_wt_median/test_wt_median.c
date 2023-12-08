#define PROG_NAME "test_wt_median"
#define PROG_DESC "test of {wt_median.h} and {wt_median_window.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-11-25 17:19:50 by stolfi */ 
/* Created on 2012-03-04 by J. Stolfi, UNICAMP */

#define test_wt_median_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <rn.h>
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

void twtm_test_median_sorted(int32_t ns_max);
  /* Tests the procedure {wt_median_sorted} with up to {ns_max} elements. */

void twtm_test_median_unsorted(int32_t n_max);
  /* Tests the procedure {wt_median_unsorted} with up to {n_max} samples and weights. */

void twtm_test_gather_samples(int32_t nx_max, int32_t nw_max);
  /* Tests the procedure {wt_median_gather_samples} with up to {nx_max}
    samples and {nw_max} weights and sample indices.
    Also tests {wt_median_index_set_sort}, {wt_median_index_set_quick_sort},
    and {wt_median_index_set_insertion_sort}. */
  
void twtm_test_index_set_update(int32_t nx_max, int32_t nw_max);
  /* Tests the procedure {wt_median_window_index_set_update} with up to {nx_max}
    samples and {nw_max} weights and sample indices. */

int32_t twtm_choose_num_samples(int32_t ktry, int32_t nx_max);
  /* Chooses a number of samples {nx} in {1..nx_max}. */

void twtm_choose_samples(int32_t nx, double x[]);
  /* Fills {x[0..nx-1]} with random sample values.  */
    
int32_t twtm_choose_window_size(int32_t ktry, int32_t nw_max, int32_t nx);
  /* Chooses a number of weights (the window size) {nw} in {1..min(nx,nw_max)}. */

void twtm_choose_weights(int32_t nw, int32_t w[]);
  /* Fills {w[0..nw-1]} with random non-negative weights. Many will be zero, but the
    sum will be positive. */

int32_t twtm_choose_indices(int32_t ktry, int32_t nx, int32_t nk, int32_t kx[]);
  /* Fills {kx[0..nk-1]} with a set of {nk} consecutive indices {ix .. ix+nk-1}
    in the range {0..nx-1}, in random order.  Returns the lowest index {ix}.
    Requires {1 <= nk <= nx}. */

void twtm_check_gather_samples
  ( int32_t nx,
    double x[],
    int32_t ix,
    int32_t nw,
    int32_t w[],
    int32_t ns,
    double xs[],
    int32_t ws[],
    int32_t kx[]
  );
  /* Checks whether the vectors {xs[0..ns-1]} and {ws[0..ns-1]} are proper condensations of 
    the samples {x[kx[0..nw-1]]}.  Assumes that {w[[kx[i]-ix]} is the weight of {x[kx[i]]}
    for {i} in {0..nw-1}, and that {kx[0..nw-1]} is a permutation of {ix..ix+nw-1]}. */
    
void twtm_check_median(int32_t ns, bool_t sorted, double xs[], int32_t ws[], bool_t interp, double xm, bool_t verbose);
  /* Checks whether {vm} is the correct output of {wt_median_sorted(ns,xs,ws,interp)}. */
  
void twtm_check_index_set_update
  ( int32_t nk0,
    int32_t kx0[],
    int32_t nx,
    int32_t ix,
    int32_t nk1,
    int32_t kx1[],
    int32_t np
  );
  /* Assumes that whether {wt_median_window_index_set_update} applied to the index set {kx0[0..nk0-1]} with 
    parameters {nx, ix, nk1} resulted in {kx1[0..nk1-1]} with {np} elements declared "new". 
    Checks if that is what was expected. */

void twtm_prxsws(int32_t ns, char *xname, double xs[], char *wname, int32_t ws[]);  
  /* Prints to stderr samples {xs[0..ns-1]} and weights
    {ws[0..ns-1]} on one line, with names {xname} and {wname}
    respectively. */

void twtm_prkxw(char *xname, int32_t nx, double x[], char *wname, int32_t nw, int32_t w[], char *kname, int32_t kx[]);
  /* Prints to stderr the indices {kx[0..nw-1]}, the samples {x[kx[0..nw-1]]}, and the weights
    {w[0..nw-1]} on one line, with names {xname} and {wname}
    respectively. */

void twtm_prxw(char *xname, int32_t nx, double x[], int32_t ix, char *wname, int32_t nw, int32_t w[]);
  /* Prints to stderr the samples {x[ix..ix+nw-1]} and weights
    {w[0..nw-1]} on one line, with names {xname} and {wname}
    respectively. */

void twtm_prkx(int32_t nk, int32_t kx[], int32_t np);
  /* Prints to stderr the perm {kx[0..nw-1]} on one line.
    Prints a '|' before {kx[np]}, if {np < nw}. */

int32_t main (int32_t argc, char **argv)
  {
    twtm_test_median_sorted(1000);
    twtm_test_gather_samples(1000, 30);
    twtm_test_median_unsorted(1000);
    twtm_test_index_set_update(1000, 30);
    return 0;
  }
  
void twtm_test_median_sorted(int32_t ns_max)
  { fprintf(stderr, "> %s ns_max = %d\n", __FUNCTION__, ns_max);
    double *xs = talloc(ns_max, double);
    int32_t *ws = talloc(ns_max, int32_t);
    for (int32_t ktry = 0; ktry < 200; ktry++)
      { bool_t verbose = (ktry < 10);
        if (verbose) { fprintf(stderr, "------------------------------------------------------------\n"); }
        
        /* Choose {ns}: */
        int32_t ns = twtm_choose_num_samples(ktry, ns_max);
        if (verbose) { fprintf(stderr, "  ns = %d\n", ns); }
        
        /* Generate {ns} distinct samples in order with nonzero weights: */
        for (int32_t ks = 0; ks < ns; ks++)
          { xs[ks] = (ks == 0 ? dabrandom(-ns/2, +ns/2) : xs[ks-1]) + exp(dabrandom(-15.0,0.0));
            ws[ks] = abrandom(1, 20);
          }
        if (verbose && (ns <= 10)) { twtm_prxsws(ns, "xs", xs, "ws", ws); }
        
        /* Compute median with and without interpolation: */
        for (int32_t ki = 0; ki < 2; ki++)
          { bool_t interp = (ki > 0);
            if (verbose) { fprintf(stderr, "  interp = %c\n", "FT"[interp]); }
            double xm = wt_median_sorted(ns, xs, ws, interp);
            if (verbose) { fprintf(stderr, "  xm = %+14.8f = %24.15e\n", xm, xm); }
            twtm_check_median(ns, TRUE, xs, ws, interp, xm, verbose);
            if (verbose) { fprintf(stderr, "\n"); }
          }
      }
    fprintf(stderr, "< %s\n", __FUNCTION__);
  }
         
void twtm_test_gather_samples(int32_t nx_max, int32_t nw_max)
  { fprintf(stderr, "> %s  nx_max = %d nw_max = %d\n", __FUNCTION__, nx_max, nw_max);

    /* Raw samples and weights: */
    double *x = talloc(nx_max, double);
    int32_t *w = talloc(nw_max, int32_t);
    
    /* Gathered samples and weights: */
    double *xs = talloc(nw_max, double);
    int32_t *ws = talloc(nw_max, int32_t);
    
    /* Index set: */
    int32_t *kx = talloc(nw_max, int32_t);

    for (int32_t ktry = 0; ktry < 200; ktry++)
      { bool_t verbose = (ktry < 10);
        if (verbose) { fprintf(stderr, "------------------------------------------------------------\n"); }

        /* Choose the samples and weights: */
        int32_t nx = twtm_choose_num_samples(ktry, nx_max);
        twtm_choose_samples(nx, x);
        int32_t nw = twtm_choose_window_size(ktry, nw_max, nx);
        twtm_choose_weights(nw, w);
        
        /* Fill the index set: */
        int32_t ix = twtm_choose_indices(ktry, nx, nw, kx);

        /* Choose sort procedure: */
        int32_t which_sort = (ktry/2) % 3;
        if (which_sort == 0)
          { int32_t np = abrandom(0, nw);
            if (verbose) { fprintf(stderr, "  sort (auto, np = %d)\n", np); }
            wt_median_index_set_sort(nx, x, nw, kx, np);
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
        int32_t ns = wt_median_gather_samples(nx, x, ix, nw, w, kx, xs, ws);
        if (verbose && (ns <= 10)) { twtm_prxsws(ns, "xs", xs, "ws", ws); }
        
        /* Check the result: */
        twtm_check_gather_samples(nx, x, ix, nw, w, ns, xs, ws, kx);
      }
    fprintf(stderr, "< %s\n", __FUNCTION__);
  }

void twtm_test_median_unsorted(int32_t n_max)
  { fprintf(stderr, "> %s ns_max = %d\n", __FUNCTION__, n_max);
    
    /* data vectors: */
    double *x = talloc(n_max, double);
    int32_t *w = talloc(n_max, int32_t);
    
    /* Work vectors: */
    double *xs = talloc(n_max, double);
    int32_t *ws = talloc(n_max, int32_t);
    int32_t *kx = talloc(n_max, int32_t);
    
    for (int32_t ktry = 0; ktry < 200; ktry++)
      { bool_t verbose = (ktry < 10);
        if (verbose) { fprintf(stderr, "------------------------------------------------------------\n"); }
        
        /* Choose number of samples {n}: */
        int32_t n = twtm_choose_num_samples(ktry, n_max);
        assert((n > 0) && (n <= n_max));
        if (verbose) { fprintf(stderr, "  n = %d\n", n); }
        
        /* Generate {n} distinct samples with non-negative weights: */
        twtm_choose_samples(n, x);
        twtm_choose_weights(n, w);
        if (verbose && (n <= 10)) { twtm_prxw("x", n, x, 0, "w", n, w); }
            
        for (int32_t ki = 0; ki < 2; ki++)
          { bool_t interp = (ki > 0);
            if (verbose) { fprintf(stderr, "  interp = %c\n", "FT"[interp]); }
            
            /* Decide whether to use external or internal tables: */
            double *xs_try = ((ktry & 9) != 0 ? xs : NULL);
            int32_t *ws_try = ((ktry & 18) != 0 ? ws : NULL);
            int32_t *kx_try = ((ktry & 36) != 0 ? kx : NULL);

            /* Compute median: */
            int32_t ns = -1;
            double xm = NAN;
            xm = wt_median_unsorted(n, x, w, interp, &ns, xs_try, ws_try, kx_try);
            if (verbose) { fprintf(stderr, "  xm = %+14.8f = %+24.15e ns = %d\n", xm, xm, ns); }
            twtm_check_median(n, FALSE, x, w, interp, xm, verbose);
            if ((xs_try == xs) && (ws_try == ws))
              { if (verbose && (ns <= 10)) { twtm_prxsws(ns, "xs", xs, "ws", ws); }
                twtm_check_median(ns, TRUE, xs, ws, interp, xm, verbose);
              }
            if (verbose) { fprintf(stderr, "\n"); }
          }
      }
    fprintf(stderr, "< %s\n", __FUNCTION__);
  }

void twtm_test_index_set_update(int32_t nx_max, int32_t nw_max)
  {
    fprintf(stderr, "> %s  nx_max = %d nw_max = %d\n", __FUNCTION__, nx_max, nw_max);

    /* Index sets: */
    int32_t *kx = talloc(nw_max, int32_t);
    int32_t *kx_save = talloc(nw_max, int32_t);

    for (int32_t ktry = 0; ktry < 200; ktry++)
      { bool_t verbose = (ktry < 10);

        if (verbose) { fprintf(stderr, "------------------------------------------------------------\n"); }

        /* Choose the samples and two windows: */
        int32_t nx = twtm_choose_num_samples(ktry, nx_max);
        int32_t nk0 = twtm_choose_window_size(ktry, nw_max, nx);
        int32_t nk1 = twtm_choose_window_size(ktry, nw_max, nx);

        if (verbose) { fprintf(stderr, "> nx = %d nk0 = %d nk1 = %d\n", nx, nk0, nk1); }
        demand((0 <= nk0) && (nk0 <= nx), "window size {nk0} too big for sample vector");
        demand((0 <= nk1) && (nk1 <= nx), "window size {nk1} too big for sample vector");
        
        /* Initial fill of the index set {kx[0..nk0-1]}, and save a copy {kx_save[1..nk0-1]}: */
        int32_t ix0 = twtm_choose_indices(ktry, nx, nk0, kx);
        if (verbose && (nk0 <= 15)) { fprintf(stderr, "  ix0 = %d\n", ix0); }
        if (verbose && (nk0 <= 15)) { twtm_prkx(nk0, kx, -1); }
        for (int32_t ik = 0; ik < nk0; ik++) { kx_save[ik] = kx[ik]; }
        
        /* Update the index set: */
        if (verbose) { fprintf(stderr, "  wt_median_window_index_set_update ...\n"); }
        int32_t ix1 = abrandom(0, nx - nk1);
        int32_t np = wt_median_window_index_set_update(nx, ix1, nk1, nk0, kx);
        if (verbose && (nk1 <= 15)) { twtm_prkx(nk1, kx, np); }

        /* Check whether it worked: */
        twtm_check_index_set_update(nk0, kx_save, nx, ix1, nk1, kx, np);
      }
  }

void twtm_check_index_set_update
  ( int32_t nk0,
    int32_t kx0[],
    int32_t nx,
    int32_t ix,
    int32_t nk1,
    int32_t kx1[],
    int32_t np
  )
  {
    /* Set {given[jx]} to true iff {jx} appears in {kx0[0..nk0-1]}: */
    bool_t given[nx];
    for (int32_t jx = 0; jx < nx; jx++) { given[jx] = FALSE; }
    for (int32_t ik = 0; ik < nk0; ik++) 
      { int32_t jxi = kx0[ik];
        demand((jxi >= 0) && (jxi < nx), "invalid index in {kx0}");
        demand(! given[jxi], "repeated indices in {kx0}");
        given[jxi] = TRUE;
      }

    /* Check {kx1[0..kn1-1]}: */
    bool_t seen[nx]; /* To check that {kx1[0..nk1-1]} are all distinct. */
    for (int32_t jx = 0; jx < nx; jx++) { seen[jx] = FALSE; }
    for (int32_t ik = 0; ik < nk1; ik++)
      { int32_t jxi = kx1[ik];
        demand((jxi >= ix) && (jxi <= ix + nk1 - 1), "{kx1[ik]} outside window index range");
        demand(! seen[jxi], "repated index in {kx1}");
        seen[jxi] = TRUE;
        if (ik < nk1-np)
          { /* The first {nk1-np} elements must come from the input set: */
            demand(given[jxi], "inconsistent {np} - new listed as old");
          }
        else
          { /* The last {np} must be new: */
            demand(! given[jxi], "inconsistent {np} - old listed as new");
          }
      }
  }

int32_t twtm_choose_num_samples(int32_t ktry, int32_t nx_max)
  { bool_t debug = FALSE;
    /* Choose {nx} and fill {x}: */
    int32_t nx = 1 + (ktry % nx_max);
    assert(nx <= nx_max);
    if (ktry > 5) { nx += abrandom(0, nx_max - nx); }
    if (debug) { fprintf(stderr, "  nx = %d\n", nx); }
    assert((nx > 0) && (nx <= nx_max));
    return nx;
  }

void twtm_choose_samples(int32_t nx, double x[])
  { 
    for (int32_t ix = 0; ix < nx; ix++) { x[ix] = dabrandom(-99,+99); }
  }

int32_t twtm_choose_window_size(int32_t ktry, int32_t nw_max, int32_t nx) 
  {
    bool_t debug = FALSE;
    /* Choose {nw} and fill {w}: */
    int32_t m = (nw_max < nx ? nw_max : nx);
    int32_t nw = 1 + (ktry % m);
    assert(nw <= m);
    if (ktry > 5) { nw += abrandom(0, m - nw); }
    if (debug) { fprintf(stderr, "  nw = %d\n", nw); }
    assert((nw > 0) && (nw <= nw_max) && (nw <= nx));
    return nw;
  }

void twtm_choose_weights(int32_t nw, int32_t w[]) 
  {
    int32_t wsum_org = 0;
    for (int32_t kw = 0; kw < nw; kw++) { w[kw] = abrandom(0,5); wsum_org += w[kw]; }
    /* Make sure that the weight sum is positive: */
    if (wsum_org == 0) { w[nw/2] = 1; wsum_org = 1; }
  }
        
int32_t twtm_choose_indices(int32_t ktry, int32_t nx, int32_t nk, int32_t kx[])
  {
    bool_t debug = FALSE;
    
    demand((nk >= 1) && (nk <= nx), "invalid {nx,nk}");
    int32_t jx0 = (nx - nk)/3;
    int32_t jx1 = jx0 + nk - 1;
    if (debug) { fprintf(stderr, "  jx0 = %d jx1 = %d\n", jx0, jx1); }
    assert(jx0 <= jx1);
    if (nk == 1)
      { kx[0] = jx0; }
    else
      { for (int32_t ik = 0; ik < nk; ik++)
          { /* Generate {kx[ik] By interpolation: */
            kx[ik] = jx0 + ik; 
            /* Shuffle it: */
            int32_t rk = abrandom(0, ik);
            int32_t t = kx[ik]; kx[ik] = kx[rk]; kx[rk] = t;
          }
      }
    return jx0;
  }

void twtm_check_median(int32_t ns, bool_t sorted, double xs[], int32_t ws[], bool_t interp, double xm, bool_t verbose)
  { 
    if (sorted)
      { /* Check strict order of {xs[0..ns-1]} and positive weights: */
        for (int32_t ks = 0; ks < ns; ks++)
          { demand(ws[ks] > 0, "invalid weights");
            if (ks > 0) { demand(xs[ks] > xs[ks-1], "samples not sorted"); }
          }
      }

    /* Compute weight sums and samples with best {F} on each side of zero: */
    bool_t is_a_sample = FALSE; /* True if {xm} is one of the {x[i]}. */
    int32_t Stot = 0; /* Sum of all weights. */
    int32_t Seq = 0;  /* Weight of {xm} if present. */
    double xlo = -INF; int32_t wlo = 0, Slo = 0;  /* largest {xs[i]} less than {xm}, its weight and below-sum. */
    double xhi = +INF; int32_t whi = 0, Shi = 0;  /* smallest {xs[i]} greter than {xm}, its weight and above-sum. */
    
    for (int32_t ks = 0; ks < ns; ks++)
      { double xsk = xs[ks];
        int32_t wsk = ws[ks];
        demand(wsk >= 0, "invalid (negative) weights");
        if (sorted) 
          { demand(wsk != 0, "invalid (zero) weights");
            if (ks > 0) { demand(xs[ks] > xs[ks-1], "samples not sorted"); }
          }
        if (wsk > 0)
          { Stot += wsk;
            if (xsk < xm)
              { Slo += wsk;
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
              { Shi += wsk;
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
               { Seq += wsk; 
                 is_a_sample = TRUE;
               }
           }
      }
    if (verbose) { fprintf(stderr, "  Stot = %d Seq = %d is_a_sample = %c\n", Stot, Seq, "FT"[is_a_sample]); }
    if (verbose) { fprintf(stderr, "  xlo = %+14.8f = %+24.15e wlo = %11d Slo = %11d\n", xlo, xlo, wlo, Slo); }
    if (verbose) { fprintf(stderr, "  xhi = %+14.8f = %+24.15e whi = %11d Shi = %11d\n", xhi, xhi, whi, Shi); }

    assert(Slo + Seq + Shi == Stot);
    assert((xlo < xm) && (xm < xhi));
    int32_t Fm = Slo - Shi;
    int32_t Flo = Fm - Seq - wlo; /* Estimated {F(xlo)}. */
    int32_t Fhi = Fm + Seq + whi; /* Estimated {F(xhi)}. */
    if (verbose) { fprintf(stderr, "  Flo = %+d  Fm = %+d  Fhi = %+d\n", Flo, Fm, Fhi); }
    if (interp)
      { if (! is_a_sample)
          { demand ((Flo < 0) && (Fhi > 0), "wrong gap for median");
            /* Check interpolation: */
            assert(Seq == 0);
            assert(isfinite(xlo) && isfinite(xhi));
            assert((Flo < 0) && (Fhi > 0));
            double f = ((double)(0 - Flo))/((double)(Fhi - Flo));
            double xm_exp = (1-f)*xlo + f*xhi;
            if (verbose) { fprintf(stderr, "  f = %14.8f = %24.15e  xm_exp = %14.8f = %24.15e\n", f, f, xm_exp, xm_exp); }
            demand(fabs(xm - xm_exp) < 1.0e-12*(fabs(xlo) + fabs(xhi)), "incorrect interpolation");
          }
        else
          { demand(Slo == Shi, "should have interpolated");
            assert(Fm == 0);
          }
      }
    else
      { demand(is_a_sample, "median should be one of the samples"); 
        demand ((fabs(Fhi) >= fabs(Fm)) && (fabs(Flo) >= fabs(Fm)), "not optimal choice for median");
      }
  }
        
void twtm_check_gather_samples
  ( int32_t nx,
    double x[],
    int32_t ix,
    int32_t nw,
    int32_t w[],
    int32_t ns,
    double xs[],
    int32_t ws[],
    int32_t kx[]
  )
  { 
    int32_t wsum_org = 0; /* Sum of all given weights. */
    for (int32_t kw = 0; kw < nw; kw++) { wsum_org += w[kw]; }
    int32_t wsum_chk = 0; /* Sum of the condensed weights. */
    double xs_prev = -INF;
    for (int32_t ks = 0; ks < ns; ks++)
      { double xsk = xs[ks];
        demand(isfinite(xsk), "condensed sample is infinite or {NAN}");
        demand(xsk > xs_prev, "condensed samples out of order");
        int32_t wsk_chk = 0;
        for (int32_t ik = 0; ik < nw; ik++)
          { int32_t jx = kx[ik];
            demand((jx >= 0) && (jx < nx), "invalid index in {kx}");
            int32_t iw = kx[ik] - ix;
            demand((iw >= 0) && (iw < nw), "index in {kx} out of window range");
            if (x[jx] == xsk) { wsk_chk += w[iw]; }
          }
        demand(wsk_chk == ws[ks], "condensed sample weight mismatch");
        demand(wsk_chk > 0, "condensed samples with zero weight");
        wsum_chk += ws[ks];
        xs_prev = xsk;
      }
    demand(wsum_chk == wsum_org, "lost weight");
  }

void twtm_prxsws(int32_t ns, char *xname, double xs[], char *wname, int32_t ws[])      
  { fprintf(stderr, "  %s[0..%d] = ", xname, ns-1);
    for (int32_t k = 0; k < ns; k++) { fprintf(stderr, " %+14.8f", xs[k]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s[0..%d] = ", wname, ns-1);
    for (int32_t k = 0; k < ns; k++) { fprintf(stderr, " %14d", ws[k]); }
    fprintf(stderr, "\n");
  }
       
void twtm_prkxw(char *xname, int32_t nx, double x[], char *wname, int32_t nw, int32_t w[], char *kname, int32_t k[])      
  { fprintf(stderr, "  nx = %d  nw = %d\n", nx, nw);
    fprintf(stderr, "  %s[0..%d] =     ", kname, nw-1);
    for (int32_t j = 0; j < nw; j++) { fprintf(stderr, " %14d", k[j]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s[%s[0..%d]] =  ", xname, kname, nw-1);
    for (int32_t j = 0; j < nw; j++) { fprintf(stderr, " %+14.8f", x[k[j]]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s[0..%d] =     ", wname, nw-1);
    for (int32_t j = 0; j < nw; j++) { fprintf(stderr, " %14d", w[j]); }
    fprintf(stderr, "\n");
  }
     
void twtm_prxw(char *xname, int32_t nx, double x[], int32_t ix, char *wname, int32_t nw, int32_t w[])      
  { int32_t jx = ix+nw-1;
    fprintf(stderr, "  n = %d\n", nw);
    fprintf(stderr, "  %s[%d..%d] = ", xname, ix, jx);
    for (int32_t k = 0; k < nw; k++) { fprintf(stderr, " %+14.8f", x[ix + k]); }
    fprintf(stderr, "\n");
    fprintf(stderr, "  %s[0..%d] = ", wname, nw-1);
    for (int32_t k = 0; k < nw; k++) { fprintf(stderr, " %14d", w[k]); }
    fprintf(stderr, "\n");
  }
        
void twtm_prkx(int32_t nk, int32_t kx[], int32_t np)
  { fprintf(stderr, "  n = %d np = %d\n", nk, np); 
    fprintf(stderr, "  kx ="); 
    for (int32_t k = 0; k < nk; k++) 
      { if (k == nk-np) { fprintf(stderr, " |"); }
        fprintf(stderr, " %d", kx[k]);
      }
    fprintf(stderr, "\n");
  }
