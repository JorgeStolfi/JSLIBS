/* Last edited on 2023-11-25 17:20:52 by stolfi */

void test_basic(int32_t nx, int32_t nw_max);
  /* Tests the correctness with {nx} elements and up to {nw_max}
    window widths. */

void test_median_two_win_sizes(int32_t nx, double x[], int32_t nw0, bool_t interp);
  /* Tests the correctness of {wt_median_index_set_update}, 
    {wt_median_index_set_sort}, }wt_median_gather_samples} {wt_median},
    with elements {x[0..nx-1]} and window widths alternating between {nw0}
    and {nw0+5}. */

void test_median_in_gap(int32_t nw, bool_t interp);
  /* Tests the correctness of {wt_median_in_window} in the special case 
    when the median falls between two elements. */


void twtm_check_median(int32_t nx, double x[], int32_t ix, int32_t nw, int32_t w[], bool_t interp, double xm, bool_t verbose);
  /* Checks whether {vm} is the correct median. */

    test_basic(1000, 30);
    /* test_timing(1000, 100); */

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
      { for (bool_t interp = FALSE; interp <= TRUE; interp++)
          { test_median_two_win_sizes(nx, x, nw0, interp); 
            test_median_in_gap(nw0, interp);
          }
      }
  }

void test_median_two_win_sizes(int32_t nx, double x[], int32_t nw0, bool_t interp)

void twtm_text_running_median()
  {

    double *wf0 = rn_alloc(nw0); wt_table_hann_fill(nw0, 0.0, wf0, NULL);
    int32_t *wi0 = talloc(nw0, int32_t);
    int32_t wi0_sum = wt_table_quantize(nw0, wf0, 100, 100000, wi0);
    demand(wi0_sum <= wt_median_WSUM_MAX, "sum of weights {wi1} too big");
    
    double *wf1 = rn_alloc(nw1); wt_table_fill_hann(nw1, 0.0, wf1, NULL);
    int32_t *wi1 = talloc(nw1, int32_t);
    int32_t wi1_sum = wt_table_quantize(nw1, wf1, 100, 100000, wi1);
    demand(wi1_sum <= wt_median_WSUM_MAX, "sum of weights {wi1} too big");

    int32_t kx[nw1];  /* Indices in window, sorted by {x} value. */
    int32_t nk = 0; /* Current indices in {kx} are {kx[0..nk-1]}. */

    bool_t seen[nw1]; /* To check whether {kx} is a permutation. */
    
    double xs[nw1];  /* Sorted and consolidated window samples. */
    int32_t ws[nw1]; /* Corresponding condensed weights. */
    
    int32_t trial = 0;
    int32_t ix = 0; /* Index of first sample in window. */
    while (ix < nx-1)
      { /* Select the window width {nw} and weight table {w}: */
        int32_t nw_max = nx - ix;
        assert(nw_max >= nw0);
        int32_t nw = ((trial % 4) < 2 ? nw0 : nw1);
        if (nw > nw_max) { nw = nw0; }
        assert((nw <= nw_max) && ((nw == nw0) || (nw == nw1)));

        int32_t jx = ix+nw-1; /* Index of last sample in window. */
        assert(jx < nx);
        
        int32_t *w = (nw == nw0 ? wi0 : wi1);
        if (verbose) 
          { fprintf(stderr, "  nx = %d nw = %d", nx, nw);
            fprintf(stderr, "  window ix = {%d..%d} interp = %c \n", ix, jx, "FT"[interp]);
          }
        if (verbose && (nw <= 9)) { twtm_prxw("x", nx, x, ix, "w", nw, w); }
          
        /* Update the index set: */
        if (verbose) { fprintf(stderr, "  wt_median_index_set_update ...\n"); }
        int32_t np = wt_median_index_set_update(nx, ix, nw, nk, kx);
        if (verbose  && (nw <= 15)) { twtm_prkx(nw, kx, np); }
        nk = nw;
        
        /* Sort the index set: */
        if (verbose) { fprintf(stderr, "  wt_median_index_set_sort ...\n"); }
        wt_median_index_set_sort(nx, x, nw, kx, np);
        if (verbose  && (nw <= 15)) { twtm_prkx(nw, kx, nw); }
        
        /* Extract and condense the samples and weights: */
        if (verbose) { fprintf(stderr, "  wt_median_gather_samples ...\n"); }
        int32_t ns = wt_median_gather_samples(nx, x, nw, w, kx, xs, ws);
        if (verbose && (nw <= 9)) { twtm_prxw("xs", ns, xs, 0, "ws", ns, ws); }
        
        /* Compute median {xm}: */
        if (verbose) { fprintf(stderr, "  wt_median ...\n"); }
        double xm = wt_median(ns, xs, ws, interp);
        fprintf(stderr, "  xm = %+10.7f\n", xm);
         
        /* Check if {xm} looks like a median: */
        twtm_check_median(nx, x, ix, nw, w, interp, xm, verbose);

        /* Check {kx[0..nw-1]}: */
        for (int32_t k = 0; k < nw; k++) { seen[k] = FALSE; }
        for (int32_t k = 0; k < nw; k++)
          { int32_t rx = kx[k];
            demand((rx >= ix) && (rx <= jx), "{kx[i]} outside window index range");
            demand(! seen[rx - ix], "repated index in {kx}");
            seen[rx - ix] = TRUE;
          }

      /* Advance {ix}: */
      int32_t ix_step = ((ix < 2*nw0) || (jx >= nx - 2*nw0));
      ix += ix_step;
      trial++;
      fprintf(stderr, "\n");
    }
  }  
void test_median_in_gap(int32_t nw, bool_t interp)
  { int32_t hw = nw/2; 
    int32_t nx = nw + 30;
    double x[nx];
    double xmin = 0.0, xmax = 300.0;
    for (int32_t ix = 0; ix < nx; ix++) { x[ix] = dabrandom(xmin, xmax); }
    int32_t w[nw];
    double xlo = 100.0,    xhi = 200.0;  /* Median will range from {xlo} to {xhi}. */
        
    int32_t nt = 10; /* Will try {nt+1} median positions. */
    for (int32_t t = 0; t <= nt; t++)
      { 
        /* Specify {w[0..nw]} symmetric about center: */
        int32_t Slo = 0; /* Sum of weights of window samples {<= xlo}. */
        int32_t Shi = 0; /* Sum of weights of window samples {>= xhi}. */
        for (int32_t k0 = 0; k0 < hw; k0++)
          { int32_t k1 = nw-1-k0;
            int32_t wk = abrandom(10,99);
            w[k0] = wk; Slo += wk;
            w[k1] = wk; Shi += wk;
          }
        /* If the window width is odd, set the extra window weight to 0: */
        if ((nw & 1) == 1) { w[hw] = 0; }

        int32_t ix = (nx - nw)/2;    /* Index of start of window. */
        int32_t jx = ix+nw-1;        /* Index of end of window. */

        /* Specify {x[ix..jx]} so that half are in {[xmin_xlo-1]}  and half in {[xhi+1_xmax]}: */
        for (int32_t k0 = 0; k0 < hw; k0++)
          { int32_t k1 = nw-1-k0;
            x[ix + k0] = dabrandom(xmin, xlo-1);
            x[ix + k1] = dabrandom(xhi+1, xmax);
          }
        /* If the window width is odd, set the extra window elem to someting in {[xlo_xhi]}. It should be ignored. */
        if ((nw & 1) == 1) { x[ix + hw] = (sqrt(2)*xlo + xhi)/(sqrt(2) + 1); }

        /* Make sure that exactly one sample is {xlo} and one is {xhi}: */
        x[ix] = xlo;
        x[jx] = xhi;
        
        /* The median should be strictly halfway between {xlo} and {xhi}. */
        assert(Slo == Shi);

        /* Perturb weights {w[0],w[1]} so as to preserve {Slo}: */
        int32_t dw0 = abrandom(1,w[0]-1);
        w[0] -= dw0; w[1] += dw0;
        assert((w[0] > 0) && (w[1] > 0));

        /* Perturb weights {w[nw-1],w[nw-2]} so as to preserve {Shi}: */
        int32_t dw1 = abrandom(1,w[nw-1]-1);
        w[nw-1] -= dw1; w[nw-2] += dw1;
        assert((w[nw-1] > 0) && (w[nw-2] > 0));

        /* Vary weights {w[1],w[nw-2]} so that the median ranges in {[xlo_xhi]}: */ 
        if (t == 0)
          { w[1] += w[0]; Slo += w[0]; }
        else if (t == nt)
          { w[nw-2] += w[nw-1]; Shi += w[nw-1]; }
        else
          { int32_t dw0 = abrandom(1, w[0]);
            int32_t dw1 = abrandom(1, w[nw-1]);
            w[1] += dw0; Slo += dw0;
            w[2] += dw1; Shi += dw1;
            /* Median should still be strictly between {xlo} and {xhi}: */
            assert(Slo == Shi);
          }
          
        int32_t Flo = (Slo - w[0]) - Shi;     /* Value of {F(xlo)}. */
        int32_t Fhi = Slo - (Shi - w[nw-1]);  /* Value of {F(xhi)}. */
        
        /* Compute expected median position: */
        double xm_exp; /* Expected median: */
        if (Flo == 0)
          { xm_exp = xlo; }
        else if (Fhi == 0)
          { xm_exp = xhi; }
        else
          { assert((Flo < 0) && (Fhi > 0));
            if (interp)
              { /* Estimate zero of {F} by ffine interpolation: */
                double f = ((double)Fhi)/((double)(Fhi - Flo)); 
                xm_exp = (1-f)*xlo + f*xhi; 
              }
            else
              { /* Choose either {xlo} or {xhi}.  Tie breaking may be wrong: */
                xm_exp = (abs(Flo) < abs(Fhi) ? xlo : xhi);
              }
          }

        fprintf(stderr, "  Slo = %12d wlo = %d Shi = %12d whi = %d\n", Slo, w[0], Shi, w[nw-1]);
        fprintf(stderr, "  F(xlo) = %+12d  F(xhi) = %+12d\n", Flo, Fhi);

        /* Randomly permute samples and their weights: */
        for (int32_t ki = 0; ki < nw; ki++)
          { int32_t kj = abrandom(ki, nw-1);
            if (ki != kj) 
              { int32_t tw = w[ki]; w[ki] = w[kj]; w[kj] = tw; 
                double  tx = x[ix + ki]; x[ix + ki] = x[ix + kj]; x[ix + kj] = tx; 
              }
          }
        
        /* Compute median: */
        int32_t ns;
        double xs[nw];
        int32_t ws[nw];
        double xm_cmp = wt_median_in_window(nx, x, ix, nw, w, interp, 0, NULL, &ns, xs, ws);
        
        /* Accept either {xlo} or {xhi} in case of tie: */
        if ((! interp) && (Flo < 0) && (Flo == -Fhi))
          { if ((xm_cmp == xlo) && (xm_exp == xhi))
              { xm_exp = xlo; }
            if ((xm_cmp == xhi) && (xm_exp == xlo))
              { xm_exp = xhi; }
          }
          
        fprintf(stderr, "  xm_exp = %24.15e xm_cmp = %24.15e\n", xm_exp, xm_cmp);
        demand(fabs(xm_exp - xm_cmp) < 1.0e-10, "wrong median");
      }
  }
        
void twtm_check_median(int32_t nx, double x[], int32_t ix, int32_t nw, int32_t w[], bool_t interp, double xm, bool_t verbose)
  { int32_t Slo = 0, Seq = 0, Shi = 0, Stot = 0;
    bool_t isanx = FALSE; /* True if {xm} is one of the {x[i]}. */
    /* Closest sample values to {xm} with nonzero weight, and their total weights: */
    double xlo = -INF; int32_t wlo = 0; /* Largest {x} value less than {xm}. */
    double xhi = +INF; int32_t whi = 0; /* Smallest {x} value greater than {xm}. */
    for (int32_t k = 0; k < nw; k++)
      { double xk = x[ix + k];
        int32_t wk = w[k];
        Stot += wk;
        if (wk != 0)
          { if (xk < xm) 
              { Slo += wk; 
                if (xk > xlo) { xlo = xk; wlo = wk; } else if (xk == xlo) { wlo += wk; }
              }
            else if (xk > xm)
              { Shi += wk; 
                if (xk < xhi) { xhi = xk; whi = wk; } else if (xk == xhi) { whi += wk; }
              }
            else 
              { Seq += wk; }
          }
        if (xm == xk) { isanx = TRUE; }
      }
    if (verbose) { fprintf(stderr, "  Slo = %11d  Seq = %11d Shi = %11d isanx = %c\n", Slo, Seq, Shi, "FT"[isanx]); }
    if (verbose) { fprintf(stderr, "  xlo = %24.15e wlo = %11d  xhi = %24.15e whi = %11d  \n", xlo, wlo, xhi, whi); }
    assert(Slo + Seq + Shi == Stot);
    assert((xlo < xm) && (xm < xhi));
    int32_t Fm = Slo - Shi;
    int32_t Fhi = Fm + Seq + whi; /* Estimated {F(xhi)}. */
    int32_t Flo = Fm - Seq - wlo; /* Estimated {F(xlo)}. */
    demand ((fabs(Fhi) >= fabs(Fm)) && (fabs(Flo) >= fabs(Fm)), "not optimal choice for median");
    if (! isanx)
      { demand(interp, "not {interp} - median should be one of the samples");
        assert(Seq == 0);
        /* !!! Should check interpolation... !!! */
      }
  }
 
