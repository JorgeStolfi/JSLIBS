#define PROG_NAME "test_wt_table"
#define PROG_DESC "test of {wt_table.h}, {wt_table_*}.h"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-20 06:56:01 by stolfi */ 
/* Created on 2012-03-04 by J. Stolfi, UNICAMP */

#define test_hermite3_COPYRIGHT \
  "Copyright © 2017  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <jsfile.h>
#include <jsprintf.h>
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

void test_basics(bool_t verbose);
  /* Basic consistency tests. */

void test_wt_table_print(uint32_t n, bool_t verbose);
  /* Generates a table with and prints it with {wt_table_print}. */

void test_wt_table_avg__wt_table_var(uint32_t n, bool_t verbose);
  /* Checks {wt_table_avg}, {wt_table_var} for a table of size {n}. */

void test_wt_table_normalize_sum__wt_table_check_normalization(uint32_t n, bool_t verbose);
  /* Checks {wt_table_normalize_sum} for a table of size {n}. */

void test_wt_table_convolution(uint32_t n, bool_t verbose);
  /* Check {wt_table_convolution} for two tables of size about {n} and 
    a random stride. */

void test_wt_table_gaussian_entry__wt_table_gaussian_loss(uint32_t n, double sigma, bool_t verbose);
  /* Checks {wt_table_gaussian_loss}  for the given {n} and {sigma}. */
    
void test_wt_wt_table_gaussian_make__table_gaussian_fill(double sigma, double maxLoss, bool_t verbose);
  /* Check {wt_table_gaussian_fill}, {wt_table_gaussian_make} for the given
    {sigma,maxLoss,norm}. */
    
void test_wt_table_kind_to_string__wt_table_make__wt_table_make_fill
  ( uint32_t n,
    wt_table_kind_t kind,
    double parm,
    bool_t norm,
    bool_t verbose
  );
  /* Check {wt_table_fill}, {wt_table_make} for the given{kind} and parameter {parm}
    and normalization option {norm}. */

void test_wt_table_quantize(uint32_t n, bool_t verbose);
  /* Tests {wt_table_quantize}. */

int32_t main (int32_t argc, char **argv)
  {
    test_basics(TRUE);
    for (int32_t n = 1; n <= 100; n = 3*n/2+1)
      { bool_t verbose = (n < 10);
        test_wt_table_print(n, verbose);
        test_wt_table_avg__wt_table_var(n, verbose);
        test_wt_table_normalize_sum__wt_table_check_normalization(n, verbose);
        test_wt_table_convolution(n, verbose);
        
        double sigma = n/5.0;
        test_wt_table_gaussian_entry__wt_table_gaussian_loss(n, sigma, verbose);
        test_wt_wt_table_gaussian_make__table_gaussian_fill(sigma, 0.0001, verbose);
        test_wt_wt_table_gaussian_make__table_gaussian_fill(sigma, 0.01, verbose);
        
        for (int32_t inorm = 0; inorm < 2; inorm++)
          { bool_t norm = (inorm > 0);

            test_wt_table_kind_to_string__wt_table_make__wt_table_make_fill(n, wt_table_kind_GAUSSIAN, 0.25*n, norm, verbose);
            test_wt_table_kind_to_string__wt_table_make__wt_table_make_fill(n, wt_table_kind_GAUSSIAN, 0.0, norm, verbose);
            test_wt_table_kind_to_string__wt_table_make__wt_table_make_fill(n, wt_table_kind_TRIANGULAR, 0.0, norm, verbose);
            test_wt_table_kind_to_string__wt_table_make__wt_table_make_fill(n, wt_table_kind_BINOMIAL, 0.0, norm, verbose);
            test_wt_table_kind_to_string__wt_table_make__wt_table_make_fill(n, wt_table_kind_UNIFORM, 17.0, norm, verbose);
            test_wt_table_kind_to_string__wt_table_make__wt_table_make_fill(n, wt_table_kind_HANN, 0.0, norm, verbose);
            test_wt_table_kind_to_string__wt_table_make__wt_table_make_fill(n, wt_table_kind_HANN, 1.0, norm, verbose);
            test_wt_table_kind_to_string__wt_table_make__wt_table_make_fill(n, wt_table_kind_HANN, 0.5, norm, verbose);
          }
          
        test_wt_table_quantize(n, verbose);
        
      }
    return 0;
  }

void test_basics(bool_t verbose)
  {
  
  }

void test_wt_table_print(uint32_t n, bool_t verbose)
  {
    if (! verbose) { return; }
    fprintf(stderr, "--- testing {wt_table_print} n = %d ---\n", n);
    double wf[n];         
    uint32_t stride;  /* Stride for partition of unit display. */
    wt_table_binomial_fill(n, wf, &stride);
    char *tname = jsprintf("binomial(%d)", n);
    wt_table_print(stderr, tname, n, wf, stride);
    free(tname);
  }


void test_wt_table_avg__wt_table_var(uint32_t n, bool_t verbose)
  {
    if (verbose){ fprintf(stderr, "--- testing {wt_table_avg}, {wt_table_var}  n = %d ---\n", n); }
    
    uint32_t nLo = n/3; 
    uint32_t nHi = n - nLo;
    if (verbose){ fprintf(stderr, "  nLo = %d  nHi = %d\n", nLo, nHi); }
    
    if ((nLo <= 0) || (nHi >= 0)) { return; }

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
    
void test_wt_table_convolution(uint32_t n, bool_t verbose)
  {
    if (verbose){ fprintf(stderr, "--- testing {wt_table_convolution}  n = %d\n", n); }
    uint32_t n1 = n;
    uint32_t n2 = uint32_abrandom(1, 2*n1);
    uint32_t stride = uint32_abrandom(1, 2*n1);
    
    if (verbose){ fprintf(stderr, "  n1 = %d n2 = %d stride = %d ---\n", n1, n2, stride); }
    
    double_vec_t wt1 = wt_table_binomial_make(n1);
    double_vec_t wt2 = wt_table_binomial_make(n2);
    double_vec_t ws = wt_table_convolution(n1, wt1.e, n2, wt2.e, stride);

    uint32_t ns = ws.ne;
    demand(ns == (n2-1)*stride + n1, "wrong size of {wt_table_convolution}");
    
    double tol = ((n1 <= 40) && (n2 <= 40) ? 0 : 1.0e-12);
    for (int32_t i = 0; i < ns; i++)
      { double sum = 0;
        for (int32_t k2 = 0; k2 < n2; k2++)
          { uint32_t k1 = i - k2*stride;
            if ((k1 >= 0) && (k1 < n1)) { sum += wt1.e[k1]*wt2.e[k2]; }
          }
        demand(fabs(ws.e[i] - sum) <= tol, "{wt_table_convolution} error");
      }
  }

void test_wt_table_normalize_sum__wt_table_check_normalization(uint32_t n, bool_t verbose)
  {
    if (verbose){ fprintf(stderr, "--- testing {wt_table_normalize_sum}, {wt_table_check_normalization} n = %d ---\n", n); }

    /* Create an uneven weight table: */
    uint32_t nLo = n/3; 
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

void test_wt_table_gaussian_entry__wt_table_gaussian_loss(uint32_t n, double sigma, bool_t verbose)
  { fprintf(stderr, "--- testing {wt_table_gaussian_entry}, {wt_table_gaussian_loss} ---\n");
    if (verbose){ fprintf(stderr, "  n = %d  sigma = %20.18f", n, sigma); }
    
    /* Compute sum of all entries in table: */
    double win = 0.0;
    for (int32_t k = 0; k < n; k++) { win += wt_table_gaussian_entry(n, k, sigma); }
    /* Compute total mass outside the table: */
    double wot = wt_table_gaussian_loss(n, sigma);
    if (verbose){ fprintf(stderr, "  inside = %20.18f  outside =  %20.18f  sum = %20.18f\n", win, wot, win+wot); }
    assert(fabs(1 - (win+wot)) < 1.0e-10);
  }
    
void test_wt_wt_table_gaussian_make__table_gaussian_fill(double sigma, double maxLoss, bool_t verbose)
  {
    if (verbose){ fprintf(stderr, "--- testing {wt_table_gaussian_fill}, {wt_table_gaussian_make}"); }
    if (verbose){ fprintf(stderr, "  n = auto  sigma = %22.18f maxLoss = %22.18f ---\n", sigma, maxLoss); }
    
    double_vec_t wm = wt_table_gaussian_make(0, sigma, maxLoss);
    uint32_t n = wm.ne;
    if (verbose){ fprintf(stderr, "  n chosen by {wt_table_gaussian_make} = %d\n", n); }
    demand(n % 2 == 1, "table should have odd length");
    double wf[n]; 
    uint32_t stride;
    wt_table_gaussian_fill(n, sigma, wf, &stride);
    if (verbose){ fprintf(stderr, "  stride returned by {wt_table_gaussian_fill} = %d\n", stride); }
    demand(stride == ((sigma == 0) || (n == 1) ? 1 : 0), "fill returned wrong {stride}");
    for (int32_t k = 0; k < n; k++) { demand(wf[k] == wm.e[k], "fill inconsistent with make"); }
    free(wm.e);
  }
     
void test_wt_table_kind_to_string__wt_table_make__wt_table_make_fill
  ( uint32_t n,
    wt_table_kind_t kind,
    double parm,
    bool_t norm,
    bool_t verbose
  )
  {
    if (verbose){ fprintf(stderr, "--- testing {wt_table_fill}, {wt_table_make}\n"); }
    
    char *tname = wt_table_kind_to_string(kind);
    if (verbose){ fprintf(stderr, "  kind = %s n = %d", tname, n); }
    if (verbose){ fprintf(stderr, "  parm = %14.8f = %24.16e norm = %c ---\n", parm, parm, "FT"[norm]); }
    
    double wf[n];   /* To be created with {fill} function, if {n} is odd. */
    double_vec_t wm = double_vec_new(0); /* To be created with {make} function. */
    uint32_t stride;  /* Stride for partition of unit check. */
    wm = wt_table_make(kind, n, parm, norm, &stride);
    uint32_t stride_fill;  /* Stride for partition of unit check. */
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
            if (verbose){ fprintf(stderr, "  fill[k] = %22.18f make[k] = %22.18f err = %24.16e\n", wf[k], wm.e[k], err); }
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
          { fprintf(stderr, "  checking partition of constant property with stride %d\n", stride);
            wt_table_check_partition_of_constant(n, wf, stride, tol_poc, die); 
          }
      }
    free(wm.e);
  }

void test_wt_table_quantize(uint32_t n, bool_t verbose)
  {
    if (verbose){ fprintf(stderr, "--- testing {wt_table_quantize} n = %d ---\n", n); }
    double wf[n];
    int32_t wi[n];
    for (int32_t pass = 0; pass < 64; pass++)
      { 
        bool_t zero_wf =     ((pass &  1) == 0); /* Normalize input weight sum to 1. */
        bool_t norm_wf =     ((pass &  2) == 0); /* Normalize input weight sum to 1. */
        bool_t flip_wf =     ((pass &  4) == 0); /* Flip sign of half the input values. */
        bool_t zero_wi_sum = ((pass &  8) == 0); /* Request output sum to be zero.*/
        bool_t huge_wi_sum = ((pass & 16) == 0); /* Request output sum to be huge.*/
        bool_t keep_nz =     ((pass & 32) == 0); /* Don't round nonzero to zero.*/

        /* Some combinations are not worth testing: */
        if (zero_wf && (flip_wf || norm_wf)) { continue; }
        if (zero_wi_sum && huge_wi_sum) { continue; }
        
        if (verbose)
          { fprintf(stderr, "  pass = %d\n", pass);
            fprintf(stderr, "  norm_wf = %c  flip_wf = %c\n", "FT"[norm_wf], "FT"[flip_wf]);
            fprintf(stderr, "  zero_wi_sum = %c  huge_wi_sum = %c \n", "FT"[zero_wi_sum], "FT"[huge_wi_sum]);
            fprintf(stderr, "  keep_nz = %c\n", "FT"[keep_nz]);
          }
        
        if (verbose){ fprintf(stderr, "  creating the input weight table ...\n"); }
    
        if (zero_wf) 
          { for (int32_t k = 0; k < n; k++) { wf[k] = 0.0; } }
        else if (n > 1) 
          { double sigma = 0.25*n; 
            wt_table_gaussian_fill(n, sigma, wf, NULL);
          }
        else
          { wt_table_binomial_fill(n, wf, NULL); }
        if (norm_wf) { wt_table_normalize_sum(n, wf); }
        if (flip_wf && (n > 1))
          { /* Flip half of the table: */
            for (int32_t k = 0; k < n; k++) 
              { double x = n*(((double)k)/((double)n-1) - 0.5);
                wf[k] *= x;
              }
          }
        if (verbose) 
          { fprintf(stderr, "  input weights:\n");
            for (int32_t k = 0; k < n; k++) 
              { fprintf(stderr, "   wf[%2u] = %+20.16f\n", k, wf[k]); }
          }
        
        uint64_t wia_sum = (zero_wi_sum ? 0 : (huge_wi_sum ? wt_table_quantize_WIA_SUM_MAX : 1000*n));
        if (verbose){ fprintf(stderr, "  quantizing the table for sum of abs weights = %lu ...\n", wia_sum); }
        uint64_t wia_sum_res = wt_table_quantize(n, wf, wia_sum, keep_nz, wi);
        if (verbose){ fprintf(stderr, "  returned sum of abs weight = %lu ...\n", wia_sum_res); }
        if (verbose) 
          { fprintf(stderr, "  output weights:\n");
            for (int32_t k = 0; k < n; k++) 
              { fprintf(stderr, "   wi[%2u] = %+20d\n", k, wi[k]); }
          }
        
        if (verbose){ fprintf(stderr, "  getting stats of the input and output weights ...\n"); }
        double wfa_sum = 0, wfa_min = +INF, wfa_max = 0;
        uint64_t wia_sum_cmp = 0;
        uint32_t wia_max = 0;
        for (int32_t k = 0; k < n; k++)
          { double wfk = wf[k];
            int32_t wik = wi[k];
            if (wfk == 0)
              { demand(wik == 0, "zero not preserved"); }
            else
              { /* Gather min, max and sum of abs input weights: */
                double wfak = fabs(wfk);
                wfa_sum += wfak;
                if (wfak < wfa_min) { wfa_min = wfak; }
                if (wfak > wfa_max) { wfa_max = wfak; }
                
                /* Gather min, max, and sum of abs output weights: */
                uint32_t wiak = (uint32_t)abs(wik);
                demand(wiak <= wt_table_quantize_WIA_MAX, "abs output weight too big");
                if (wiak == 0) 
                  { demand(! keep_nz, "{keep_nz} not honored"); }
                else
                  { if (! keep_nz) { demand(wia_sum > 0, "should have zeroed all"); }
                    demand((wfk < 0) == (wik < 0), "sign not preserved"); 
                    wia_sum_cmp += (uint32_t)wiak;
                    if (wiak > wia_max) { wia_max = wiak; }
                  }
              }
          }
        if (verbose)
          { fprintf(stderr, "  input abs weights:");
            fprintf(stderr, "  sum = %24.16e", wfa_sum);
            fprintf(stderr, "  min (nonzero) = %24.16e  max = %24.16e\n", wfa_min, wfa_max);
            fprintf(stderr, "  output abs weights:");
            fprintf(stderr, "  sum = %24lu  max = %24u\n", wia_sum, wia_max);
          }
        demand(wia_sum_res == wia_sum_cmp, "returned sum is incorrect");

        if (verbose){ fprintf(stderr, "  checking the quantization ...\n"); }
        if (zero_wi_sum) 
          { assert(wia_sum == 0);
            if (keep_nz && (wfa_max > 0)) 
              { demand(wia_sum_res != 0, "some inputs are nonzero but returned sum is zero"); }
            else
              { demand(wia_sum_res == 0, "input or requested sum of abs is zero, but returned sum of abs is not zero"); }
          }
        else if (wfa_sum == 0)
          { assert(wfa_max == 0);
            demand(wia_sum_res == 0, "input sum of abs is zero but returned sum of abs is not zero");
          }
        else
          { /* Check if scaling is roughly OK: */
            assert(wia_sum > 0);
            assert(wfa_sum > 0);
            assert(wia_sum_cmp > 0);
            assert(wia_max > 0);
            double scale = ((double)wia_sum_cmp)/wfa_sum; /* Approx scale effectively used. */
            for (int32_t k = 0; k < n; k++)
              { double wfk = wf[k];
                int32_t wik_cmp = wi[k];
                int32_t wiak_cmp = abs(wik_cmp);
                int32_t wiak_exp = (int32_t)floor(fabs(wfk)*scale + 0.5);
                if ((wfk != 0) && (wiak_exp == 0) && keep_nz) { wiak_exp = 1; }
                /* Should be scaled and rounded: */
                demand(abs(wiak_cmp - wiak_exp) <= 2, "roundoff error too big");
              }
          }
      }
  }
