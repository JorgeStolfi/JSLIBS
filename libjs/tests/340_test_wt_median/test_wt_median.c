#define PROG_NAME "test_wt_median"
#define PROG_DESC "test of {wt_median.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-11-04 00:06:11 by stolfi */ 
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
#include <wt_table.h>

#include <wt_median.h>

int32_t main(int32_t argn, char **argv);

void test_basic(int32_t nx, int32_t nw_max);
  /* Tests the correctness with {nx} elements and up to {nw_max}
    window widths. */

void test_basic_two_win_sizes(int32_t nx, double x[], int32_t nw0);
  /* Tests the correctness with elements {x[0..nx-1]} and 
    window widths varying between {nw0} (which must be odd)
    and {nw0+4}. */

void do_test(int32_t nx, double x[], int32_t nw, double w[]);
  /* Tests {wt_median()} on {x[0..nx-1]} with weights {w[0..nw-1]}. */

int32_t main (int32_t argc, char **argv)
  {
    test_basic(1000, 30);
    /* test_timing(1000, 100); */
    return 0;
  }
  
void test_basic(int32_t nx, int32_t nw_max)
  { 
    /* Create a vector {x} with a Brownian-like profile: */
    double *x = rn_alloc(nx);
    rn_throw_cube(nx, x);
    double ema = x[0];
    for (int32_t ix = 1; ix < nx; ix++)
      { ema = 0.8*ema + 0.2*x[ix];
        x[ix] = ema;
      }
    
    for (int32_t nw0 = 1; nw0 < nw_max; nw0 = 2*nw0 + 1)
      { test_basic_two_win_sizes(nx, x, nw0); }

  }

void test_basic_two_win_sizes(int32_t nx, double x[], int32_t nw0)
  { 
    /* The two window sizes: */
    demand((nw0 % 2) == 1, "nw0 must be odd");
    int32_t nw1 = nw0 + 4;
    int32_t hw0 = nw0/2;
    
    fprintf(stderr, "> %s : nx = %d nw0 = %d nw1 = %d\n", __FUNCTION__, nx, nw0, nw1);
        
    double *w0 = rn_alloc(nw0); wt_table_fill_hann(nw0, w0, TRUE);
    double *w1 = rn_alloc(nw1); wt_table_fill_hann(nw1, w1, TRUE);
    
    int32_t kx[nw1];
    bool_t seen[nw1]; /* To check whether {kx} is a permutation. */
    int32_t nk = 0;
    int32_t trial = 0;
    int32_t ix_min = hw0, ix_max = nx - 1 - hw0;
    int32_t ix = ix_min;
    while (ix < ix_max)
      { /* Select the window width {nw} and weight table {w}: */
        int32_t nw_max = 2*(int32_t)imin(ix, nx-1-ix) + 1;
        assert(nw_max >= nw0);
        int32_t nw = ((trial % 4) < 2 ? nw0 : nw1);
        if (nw > nw_max) { nw = nw0; }
        assert((nw <= nw_max) && ((nw == nw0) || (nw == nw1)));
        int32_t hw = nw/2;
        double *w = (nw == nw0 ? w0 : w1);
        bool_t interp = ((trial % 5) < 2);
        fprintf(stderr, "  nx = %d ix = %d nw = %d interp = %c \n", nx, ix, nw, "FT"[interp]);
        
        if (nw <= 9)
          { fprintf(stderr, "  x[%d..%d] = ", ix-hw, ix+hw);
            for (int32_t k = -hw; k <= +hw; k++) { fprintf(stderr, " %14.8f", x[ix+k]); }
            fprintf(stderr, "\n");
          }
          
        /* Compute median {xm}: */
        int32_t nk0 = nk;
        double xm = wt_median(nx, x, ix, nw, w, interp, &nk, kx);
        fprintf(stderr, "  xm = %+10.7f  nk = %d -> %d\n", xm, nk0, nk);
         
        /* Check if {xm} looks like a median: */
        double Slo = 0, Shi = 0;
        for (int32_t k = -hw; k < +hw; k++)
          { double xk = x[ix + k];
            if (xk < xm) { Slo += w[hw + k]; }
            if (xk > xm) { Shi += w[hw + k]; }
          }
        fprintf(stderr, "  Slo = %18.16f  Shi = %18.16f\n", Slo, Shi);
        demand((Slo < 0.5 + 1e-10) && (Shi < 0.5 + 1e-10), "not median");

        /* Check {kx[0..nw-1]}: */
        assert(nk == nw);
        for (int32_t ik = 0; ik < nw; ik++) { seen[ik] = FALSE; }
        for (int32_t ik = 0; ik < nw; ik++)
          { int32_t jx = kx[ik];
            int32_t k = jx-ix;
            demand((k >= -hw) && (k <= + hw), "{kx[i]} outside window index range");
            demand(!seen[hw + k], "repated index in {kx}");
            seen[hw + k] = TRUE;
          }
          
        /* Advance {ix}: */
        ix +=  ((ix < ix_min + 2*nw0) || (ix > ix_max - 2*nw0) ? 1 : 10); 
        trial++;
        fprintf(stderr, "\n");
      }
    
    fprintf(stderr, "< %s\n", __FUNCTION__);
  }
