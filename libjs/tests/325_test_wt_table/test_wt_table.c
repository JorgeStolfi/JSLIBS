#define PROG_NAME "test_wt_table"
#define PROG_DESC "test of {wt_table.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-11-26 05:52:40 by stolfi */ 
/* Created on 2012-03-04 by J. Stolfi, UNICAMP */

#define test_hermite3_COPYRIGHT \
  "Copyright © 2017  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <jsfile.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <affirm.h>

#include <wt_table.h>
#include <wt_table_generic.h>
#include <wt_table_gaussian.h>
#include <wt_table_hann.h>
#include <wt_table_binomial.h>
#include <wt_table_triangular.h>
#include <wt_table_uniform.h>
#include <wt_table_quantize.h>

int32_t main(int32_t argn, char **argv);

void do_test_basics(void);
  /* Basic consistency tests. */

void do_test_quantize(int32_t n);
  /* Tests {wt_table_quantize}. */

void do_test_avg_var(int32_t n);
  /* Checks {wt_table_avg}, {wt_table_var} for a table of size {n}. */

void do_test_convolution(int32_t n);
  /* Check {wt_table_convolution} for two tables of size about {n} and 
    a random stride. */

void do_test_normalize_sum(int32_t n);
  /* Checks {wt_table_normalize_sum} for a table of size {n}. */

void do_test_gaussian_loss(int32_t n, double sigma);
  /* Checks {wt_table_gaussian_loss}  for the given {n} and {sigma}. */
    
void do_test_gaussian_fill_make(double sigma, double maxLoss);
  /* Check {wt_table_gaussian_fill}, {wt_table_gaussian_make} for the given
    {sigma,maxLoss,norm}. */
    
void do_test_make_fill(int32_t n, wt_table_kind_t kind, double parm, bool_t norm);
  /* Check {wt_table_fill}, {wt_table_make} for the given{kind} and parameter {parm}
    and normalization option {norm}. */

void do_test_print(int32_t n);
  /* Generates a table with and prints it with {wt_table_print}. */

int32_t main (int32_t argc, char **argv)
  {
    do_test_basics();
    for (int32_t n = 1; n <= 100; n = 3*n/2+1)
      { do_test_print(n);
        do_test_avg_var(n);
        do_test_normalize_sum(n);
        do_test_convolution(n);
        double sigma = n/5.0;

        do_test_gaussian_loss(n, sigma);
        do_test_gaussian_fill_make(sigma, 0.0001);
        do_test_gaussian_fill_make(sigma, 0.01);
        for (int32_t inorm = 0; inorm < 2; inorm++)
          { bool_t norm = (inorm > 0);

            do_test_make_fill(n, wt_table_kind_GAUSSIAN, 0.25*n, norm);
            do_test_make_fill(n, wt_table_kind_GAUSSIAN, 0.0, norm);
            do_test_make_fill(n, wt_table_kind_TRIANGULAR, 0.0, norm);
            do_test_make_fill(n, wt_table_kind_BINOMIAL, 0.0, norm);
            do_test_make_fill(n, wt_table_kind_UNIFORM, 17.0, norm);
            do_test_make_fill(n, wt_table_kind_HANN, 0.0, norm);
            do_test_make_fill(n, wt_table_kind_HANN, 1.0, norm);
            do_test_make_fill(n, wt_table_kind_HANN, 0.5, norm);
          }
          
        do_test_quantize(n);
        
      }
    return 0;
  }

void do_test_basics(void)
  {
  
  }

void do_test_quantize(int32_t n)
  {
    bool_t verbose = (n < 10);
    fprintf(stderr, "--- testing {wt_table_quantize} n = %d ---\n", n);
    /* Create an weight table: */
    double wf[n];
    int32_t wi[n];
    for (int32_t pass = 0; pass < 16; pass++)
      { bool_t norm = ((pass & 1) == 0);
        if (n > 1) 
          { double sigma = 0.25*n; 
            wt_table_gaussian_fill(n, sigma, wf, NULL);
          }
        else
          { wt_table_binomial_fill(n, wf, NULL); }
        if (norm) { wt_table_normalize_sum(n, wf); }
        if ((n > 1) && ((pass & 2) == 0))
          { /* Flip half of the table: */
            for (int32_t k = 0; k < n; k++) 
              { double x = n*(((double)k)/((double)n-1) - 0.5);
                wf[k] *= x;
              }
          }
        if (verbose) 
          { fprintf(stderr, "  wf =");
            for (int32_t k = 0; k < n; k++) { fprintf(stderr, " %+11.8f", wf[k]); }
            fprintf(stderr, "\n");
          }
        int32_t wi_min = ((pass & 4) == 0 ? 100 : 0);
        int32_t wi_sum = ((pass & 8) == 0 ? 1000*n : (wi_min == 0 ? 0 : wt_table_quantize_WI_SUM_MAX));
        int32_t wi_sum_res = wt_table_quantize(n, wf, wi_min, wi_sum, wi);
        demand(wi_sum_res > 0, "returned sum is zero");
        /* Gather some data about {wf,wi}: */
        double wf_sum = 0, wf_min = +INF, wf_max = -INF;
        int32_t wi_sum_cmp = 0, wi_max_cmp = 0;
        for (int32_t k = 0; k < n; k++)
          { double wfk = wf[k];
            int32_t wik = wi[k];
            if (wfk == 0)
              { demand(wik == 0, "zero not preserved"); }
            else
              { double wfak = fabs(wfk);
                wf_sum += wfak;
                if (wfak < wf_min) { wf_min = wfak; }
                if (wfak > wf_max) { wf_max = wfak; }
                int32_t wiak = abs(wik);
                if (wi_min != 0) { demand(wiak >= wi_min, "{wi_min} not honored"); }
                if (wiak != 0) { demand((wfk < 0) == (wik < 0), "sign not preserved"); }
                wi_sum_cmp += wiak;
                if (wiak > wi_max_cmp) { wi_max_cmp = wiak; }
              }
          }
        demand(wi_sum_res == wi_sum_cmp, "returned sum is incorrect");
        assert(wi_max_cmp > 0);
        if (wi_max_cmp > 2*(wi_min + 1))
          { /* Check if scaling is roughly OK: */
            double scale = ((double)wi_max_cmp)/wf_max;
            for (int32_t k = 0; k < n; k++)
              { double wfk = wf[k];
                int32_t wik = wi[k];
                int32_t wiak = abs(wik);
                if ((wiak > 0) && (wiak > wi_min))
                  { /* Should be scaled and rounded: */
                    int32_t wiak_cmp = (int32_t)floor(fabs(wfk)*scale + 0.5);
                    demand(abs(wiak_cmp - wiak) <= 2, "rounding too big");
                  }
              }
          }
      }
        
    if (verbose) { fprintf(stderr, "--- end testing {wt_table_quantize} ---\n"); }
  }

void do_test_avg_var(int32_t n)
  {
    int32_t nLo = n/3; 
    int32_t nHi = n - nLo;
    if ((nLo <= 0) || (nHi >= 0)) { return; }
    
    fprintf(stderr, "--- testing {wt_table_avg,wt_table_var} n = %d ---\n", n);
    
    /* Create an uneven weight table: */
    double wf[n];
    double vLo = 3.14;
    double vHi = 2.18;
    for (int32_t k = 0; k < n; k++) { wf[k] = (k < nLo ? vLo : vHi); }
    double totLo = vLo*nLo;
    double totHi = vHi*nHi; 
    
    /* Check {wt_table_avg}: */
    double avgLo = ((double)nLo-1)/2;
    double avgHi = nLo + ((double)nHi-1)/2;
    double avgExp = (totLo*avgLo + totHi*avgHi)/(totLo + totHi);
    double avgCmp = wt_table_avg(n, wf);
    demand(! isnan(avgCmp), "{wt_table_avg} returned {NAN}");
    demand(fabs(avgCmp - avgExp) < 1.0e-10, "{wt_table_avg} failed");
    
    /* Check {wt_table_var}: */
    double varLo = (nLo-1)*(nLo+1)/12 + (avgLo - avgExp)*(avgLo - avgExp);
    double varHi = (nHi-1)*(nHi+1)/12 + (avgHi - avgExp)*(avgHi - avgExp);
    double varExp = (totLo*varLo + totHi*varHi)/(totLo + totHi);
    double varCmp = wt_table_var(n, wf, avgExp);
    demand(! isnan(varCmp), "{wt_table_var} returned {NAN}");
    demand(fabs(varCmp - varExp) < 1.0e-10, "{wt_table_var} failed");
  }
    
void do_test_convolution(int32_t n)
  {
    int32_t n1 = n;
    int32_t n2 = int32_abrandom(1, 2*n1);
    int32_t stride = int32_abrandom(1, 2*n1);
    
    fprintf(stderr, "--- testing {wt_table_convolution}");
    fprintf(stderr, " n1 = %d n2 = %d stride = %d ---\n", n1, n2, stride);
    
    double_vec_t wt1 = wt_table_binomial_make(n1);
    double_vec_t wt2 = wt_table_binomial_make(n2);
    double_vec_t ws = wt_table_convolution(n1, wt1.e, n2, wt2.e, stride);

    int32_t ns = ws.ne;
    demand(ns == (n2-1)*stride + n1, "wrong size of {wt_table_convolution}");
    
    double tol = ((n1 <= 40) && (n2 <= 40) ? 0 : 1.0e-12);
    for (int32_t i = 0; i < ns; i++)
      { double sum = 0;
        for (int32_t k2 = 0; k2 < n2; k2++)
          { int32_t k1 = i - k2*stride;
            if ((k1 >= 0) && (k1 < n1)) { sum += wt1.e[k1]*wt2.e[k2]; }
          }
        demand(fabs(ws.e[i] - sum) <= tol, "{wt_table_convolution} error");
      }
  }

void do_test_normalize_sum(int32_t n)
  {
    fprintf(stderr, "--- testing {wt_table_normalize_sum}, {wt_table_check_normalization} n = %d ---\n", n);

    /* Create an uneven weight table: */
    int32_t nLo = n/3; 
    double wf[n];
    double vLo = 3.14;
    double vHi = 2.18;
    for (int32_t k = 0; k < n; k++) { wf[k] = (k < nLo ? vLo : vHi); }
    
    /* Check {wt_table_normalize_sum}: */
    wt_table_normalize_sum(n, wf);
    double sum = 0;
    for (int32_t k = 0; k < n; k++) { sum += wf[k]; }
    demand(! isnan(sum), "{wt_table_normalize_sum} created {NAN}");
    double tol = 1.0e-12;
    demand(fabs(sum - 1.0) <= tol, "{wt_table_normalize_sum} failed");
    bool_t die = TRUE;
    wt_table_check_normalization(n, wf, tol, die);
  }

void do_test_gaussian_loss(int32_t n, double sigma)
  { fprintf(stderr, "--- testing {wt_table_gaussian_loss} n = %d  sigma = %22.18f ---\n", n, sigma);
    
    /* Compute sum of all entries in table: */
    double win = 0.0;
    for (int32_t k = 0; k < n; k++) { win += wt_table_gaussian_entry(n, k, sigma); }
    /* Compute total mass outside the table: */
    double wot = wt_table_gaussian_loss(n, sigma);
    fprintf(stderr, "gaussian sigma =  %20.18f n = %d", sigma, n);
    fprintf(stderr, "  inside = %20.18f  outside =  %20.18f  sum = %20.18f\n", win, wot, win+wot);
    assert(fabs(1 - (win+wot)) < 1.0e-10);
  }
    
void do_test_gaussian_fill_make(double sigma, double maxLoss)
  {
    fprintf(stderr, "--- testing {wt_table_gaussian_fill}, {wt_table_gaussian_make} with auto {n}");
    fprintf(stderr, " sigma = %22.18f maxLoss = %22.18f ---\n", sigma, maxLoss);
    
    double_vec_t wm = wt_table_gaussian_make(0, sigma, maxLoss);
    int32_t n = wm.ne;
    fprintf(stderr, "n = %d\n", n);
    demand(n % 2 == 1, "table should have odd length");
    double wf[n]; 
    int32_t stride;
    wt_table_gaussian_fill(n, sigma, wf, &stride);
    demand(stride == ((sigma == 0) || (n == 1) ? 1 : 0), "fill returned wrong {stride}");
    for (int32_t k = 0; k < n; k++) { demand(wf[k] == wm.e[k], "fill inconsistent with make"); }
    free(wm.e);
  }
     
void do_test_make_fill(int32_t n, wt_table_kind_t kind, double parm, bool_t norm)
  {
    char *tname = wt_table_kind_to_string(kind);
    fprintf(stderr, "--- testing {wt_table_fill}, {wt_table_make}\n");
    fprintf(stderr, " kind = %s n = %d", tname, n);
    fprintf(stderr, "  parm = %14.8f = %24.16e norm = %c ---\n", parm, parm, "FT"[norm]);
    
    double wf[n];   /* To be created with {fill} function, if {n} is odd. */
    double_vec_t wm = double_vec_new(0); /* To be created with {make} function. */
    int32_t stride;  /* Stride for partition of unit check. */
    wm = wt_table_make(kind, n, parm, norm, &stride);
    int32_t stride_fill;  /* Stride for partition of unit check. */
    wt_table_fill(kind, n, parm, wf, norm, &stride_fill);
    demand(stride == stride_fill, "mismatched make and fill {stride}"); 
    if (n < 20) { wt_table_print(stderr, tname, n, wf, stride_fill); }

    if (wm.ne != n)
      { fprintf(stderr, "** {wt_make_%s} returned wrong size = %d\n", tname, wm.ne); 
        assert(FALSE);
      } 
    for (int32_t k = 0; k < n; k++)
      { if (wf[k] != wm.e[k])
          { fprintf(stderr, "** discrepancy between {wt_make_%s} and {wt_fill_%s}\n", tname, tname);
            double err = wf[k] - wm.e[k];
            fprintf(stderr, "  fill[k] = %22.18f make[k] = %22.18f err = %24.16e\n", wf[k], wm.e[k], err);
            assert(FALSE);
          } 
      }
      
    if (norm)
      { double tol_norm = 1.0e-12;
        bool_t die = TRUE;
        wt_table_check_normalization(n, wf, tol_norm, die);
      }

    if (stride != 0)
      { double tol_poc = 1.0e-10;  /* Tolerance for checking values. */
        if (! norm)
          { /* Some cases are exact: */
            switch(kind)
              { case wt_table_kind_GAUSSIAN:
                  if (parm == 0) { /* Should be exact: */ tol_poc = 0.0; }
                  break;
                case wt_table_kind_BINOMIAL: 
                  if (n <= 40) 
                    { /* Should be exact, even with {norm}: */ tol_poc = 0.0; }
                  else
                    { tol_poc = 1.0e-12*wf[n/2]; }
                  break;
                case wt_table_kind_UNIFORM: 
                case wt_table_kind_TRIANGULAR: 
                  /* Exact if not normalized: */
                  tol_poc = 0.0;
                  break; 
                case wt_table_kind_HANN: 
                  if (parm == 1.0) { /* Should be exact 1.0: */ tol_poc = 0.0; }
                  break;
                default:
                  assert(FALSE);
              } 
          }
        bool_t die = TRUE; /* Abort on error. */
        if (stride_fill != 0)
          { fprintf(stderr, "checking partition of constant property with stride %d\n", stride);
            wt_table_check_partition_of_constant(n, wf, stride, tol_poc, die); 
          }
      }
    free(wm.e);
  }

void do_test_print(int32_t n)
  {
    fprintf(stderr, "--- testing {wt_table_print} n = %d ---\n", n);
    
    double wf[n];         
    int32_t stride;  /* Stride for partition of unit display. */
    wt_table_binomial_fill(n, wf, &stride);
    char *tname = NULL;
    asprintf(&tname, "binomial(%d)", n);
    wt_table_print(stderr, tname, n, wf, stride);
    free(tname);
  }

