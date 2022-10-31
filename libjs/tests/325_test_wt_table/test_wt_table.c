#define PROG_NAME "test_wt_table"
#define PROG_DESC "test of {wt_table.h}"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-30 23:19:43 by stolfi */ 
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

#include <wt_table.h>

#include <bool.h>
#include <jsfile.h>
#include <jsmath.h>
#include <affirm.h>

int32_t main(int32_t argn, char **argv);

void do_test_basics(void);
  /* Basic consistency tests. */

void do_test_avg_var(int32_t n);
  /* Checks {wt_table_avg}, {wt_table_var}. */

void do_test_normalize_sum(int32_t n);
  /* Checks {wt_table_normalize_sum}. */

void do_test_gaussian_loss(int32_t n, double sigma);
  /* Checks {wt_table_gaussian_loss}. */
    
void do_test_make_fill_gaussian(double sigma, double maxLoss, bool_t norm);
  /* Check {wt_table_fill_gaussian}, {wt_table_make_gaussian}. */
    
void do_test_make_fill(int32_t n, char *tname, bool_t norm);
  /* Check {wt_table_fill_{tname}}, {wt_table_make_{tname}}
    where {tname} is "binomial", "triangular",
    or "hann". */

void do_test_print(int32_t n, char *tname, bool_t norm);
  /* Generates aa table with {wt_table_make_{tname}(n,...,norm)} and prints it
    with {wt_table_print}. */

int32_t main (int32_t argc, char **argv)
  {
    do_test_basics();
    for (int32_t n = 1; n <= 13; n = 3*n/2+1)
      { do_test_avg_var(n);
        do_test_normalize_sum(n);
        double sigma = n/5.0;
        do_test_gaussian_loss(n, sigma);
        for (int32_t inorm = 0; inorm < 2; inorm++)
          { bool_t norm = (inorm == 1);
            do_test_print(n, "gaussian", norm);
            do_test_print(n, "binomial", norm);
            do_test_print(n, "triangular", norm);
            do_test_print(n, "hann", norm);
            do_test_make_fill_gaussian(sigma, 0.01, norm);
            do_test_make_fill_gaussian(sigma, 0.0001, norm);
            do_test_make_fill(n, "binomial", norm);
            do_test_make_fill(n, "triangular", norm);
            do_test_make_fill(n, "hann", norm);
          }
      }
    return 0;
  }

void do_test_basics(void)
  {
  
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
    
void do_test_make_fill_gaussian(double sigma, double maxLoss, bool_t norm)
  {
    fprintf(stderr, "--- testing {wt_table_fill_gaussian}, {wt_table_make_gaussian}");
    fprintf(stderr, " sigma = %22.18f maxLoss = %22.18f norm = %c---\n", sigma, maxLoss, "FT"[norm]);
    
    double_vec_t wm = wt_table_make_gaussian(sigma, maxLoss, norm);
    int32_t n = wm.ne;
    fprintf(stderr, "n = %d\n", n);
    demand(n % 2 == 1, "table shoudl have odd length");
    double wf[n]; 
    wt_table_fill_gaussian(sigma, n, wf, norm); 
    for (int32_t k = 0; k < n; k++) { assert(wf[k] == wm.e[k]); }

    if (norm)
      { bool_t die = TRUE;
        double tol = 1.0e-12;
        wt_table_check_normalization(n, wf, tol, die);
      }
    free(wm.e);
  }
     
void do_test_make_fill(int32_t n, char *tname, bool_t norm)
  {
    fprintf(stderr, "--- testing {wt_table_fill_%s}, {wt_table_make_%s}", tname, tname);
    fprintf(stderr, " n = %d norm = %c ---\n", n, "FT"[norm]);
    
    double_vec_t wm = double_vec_new(0); /* To be created with {make} function. */
    double wf[n];   /* To be created with {fill} function, if {n} is odd. */
    int32_t stride;  /* Stride for partition of unit check. */
    double tol = 1.0e-10;  /* Tolerance for unit sum and partition of unit checks. */
    
    if (strcmp(tname, "binomial") == 0)
      { wm = wt_table_make_binomial(n, norm);
        wt_table_fill_binomial(n, wf, norm);
        stride = (n == 1 ? 1 : 2); 
        tol = (n <= 40 ? 0 : 1.0e-10);  /* Shoud be exact for {n} up to 40 at least. */
      }
    else if (strcmp(tname, "triangular") == 0)
      { wm = wt_table_make_triangular(n, norm);
        wt_table_fill_triangular(n, wf, norm);
        stride = (n % 2 == 1? (n+1)/2 : 0);
      }
    else if (strcmp(tname, "hann") == 0)
      { wm = wt_table_make_hann(n, norm);
        wt_table_fill_hann(n, wf, norm); 
        stride = (n % 2 == 1? (n+1)/2 : 0);
      }
    else 
      { assert(FALSE); }

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

    if (tol < +INF)
      { bool_t die = TRUE; /* Abort on error. */
        if (norm) { wt_table_check_normalization(n, wf, tol, die); }
        /* Check partition of unit property with shift {stride}: */
        if (norm && (stride != 0))
          { fprintf(stderr, "checking partition of unity property with stride %d\n", stride);
            wt_table_check_partition_of_unity(n, wf, stride, tol, die); 
          }
      }
    free(wm.e);
  }

void do_test_print(int32_t n, char *tname, bool_t norm)
  {
    fprintf(stderr, "--- testing {wt_table_print} tname = %s n = %d norm = %c ---\n", tname, n, "FT"[norm]);
    
    double wf[n];         
    int32_t stride;  /* Stride for partition of unit display. */
    
    if (strcmp(tname, "gaussian") == 0)
      { double sigma = n/5.0;
        wt_table_fill_gaussian(sigma, n, wf, norm); 
        stride = 0; /* No partition of unit check. */
      }
    else if (strcmp(tname, "binomial") == 0)
      { wt_table_fill_binomial(n, wf, norm); 
        stride = (n == 1 ? 1 : 2);
      }
    else if (strcmp(tname, "triangular") == 0)
      { wt_table_fill_triangular(n, wf, norm);
        stride = (n % 2 == 1? (n+1)/2 : 0);
      }
    else if (strcmp(tname, "hann") == 0)
      { wt_table_fill_hann(n, wf, norm); 
        stride = (n % 2 == 1? (n+1)/2 : 0);
      }
    else 
      { assert(FALSE); }

    wt_table_print(stderr, tname, n, wf, stride);
  }

