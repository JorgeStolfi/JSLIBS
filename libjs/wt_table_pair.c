/* See wt_table_pair.h */
/* Last edited on 2024-11-20 06:51:21 by stolfi */

#define wt_table_pair_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <vec.h>
#include <bool.h>
#include <affirm.h>
#include <jsprintf.h>
#include <wt_table.h>
#include <wt_table_gaussian.h>
#include <wt_table_binomial.h>
#include <wt_table_triangular.h>

#include <wt_table_pair.h>

void wt_table_pair_make_gaussian
  ( double var0,
    double_vec_t *wt0_P, 
    char **wname0_P, 
    double var1,
    double_vec_t *wt1_P, 
    char **wname1_P,
    bool_t verbose
  )   
  { double maxLoss = 1.0e-3;
    double sigma0 = sqrt(var0);
    /* Create table with the length that it takes for the given {sigma}: */
    double_vec_t wt0 = wt_table_gaussian_make(0, sigma0, maxLoss);
    wt_table_normalize_sum(wt0.ne, wt0.e);
    uint32_t n0 = wt0.ne;
    char *wname0 = jsprintf("gaussian(n=%d,sigma=%8.6f)", n0, sigma0);
    if (verbose) { wt_table_print(stderr, wname0, wt0.ne, wt0.e, 0); }
    
    double sigma1 = sqrt(var1);
    double_vec_t wt1 = wt_table_gaussian_make(0, sigma1, maxLoss);
    wt_table_normalize_sum(wt1.ne, wt1.e);
    uint32_t n1 = wt1.ne;
    char *wname1 = jsprintf("gaussian(n=%d,sigma=%8.6f)", n1, sigma1);
    if (verbose) { wt_table_print(stderr, wname1, wt1.ne, wt1.e, 0); }
    
    (*wt0_P) = wt0; (*wname0_P) = wname0;
    (*wt1_P) = wt1; (*wname1_P) = wname1;
  }
   
void wt_table_pair_make_binomial
  ( double var0,
    double_vec_t *wt0_P, 
    char **wname0_P, 
    double var1,
    double_vec_t *wt1_P, 
    char **wname1_P,
    bool_t verbose
  )   
  { /* Variance of the distribution {w[k] = choose(n,k)} is {n/4 = r/2}, so: */
    uint32_t r0 = (uint32_t)ceil(2*var0);
    uint32_t n0 = 2*r0 + 1;
    double_vec_t wt0 = wt_table_binomial_make(n0);
    wt_table_normalize_sum(wt0.ne, wt0.e);
    char *wname0 = jsprintf("binomial(n=%d)", n0);
    if (verbose) { wt_table_print(stderr, wname0, wt0.ne, wt0.e, 0); }
    
    uint32_t r1 = (uint32_t)ceil(2*var1);
    uint32_t n1 = 2*r1 + 1;
    double_vec_t wt1 = wt_table_binomial_make(n1);
    wt_table_normalize_sum(wt1.ne, wt1.e);
    char *wname1 = jsprintf("binomial(n=%d)", n1);
    if (verbose) { wt_table_print(stderr, wname1, wt1.ne, wt1.e, 0); }
    
    (*wt0_P) = wt0; (*wname0_P) = wname0;
    (*wt1_P) = wt1; (*wname1_P) = wname1;
  }
   
void wt_table_pair_make_triangular
  ( double var0,
    double_vec_t *wt0_P, 
    char **wname0_P, 
    double var1,
    double_vec_t *wt1_P, 
    char **wname1_P,
    bool_t verbose
  )   
  { /* Variance of the distribution {w[k] = r + 1 - |r-k|} is {~r^2/4}, so: */
    uint32_t r0 = (uint32_t)ceil(sqrt(4*var0));
    uint32_t n0 = 2*r0 + 1;
    double_vec_t wt0 = wt_table_triangular_make(n0);
    char *wname0 = jsprintf("triangular(n=%d)", n0);
    wt_table_normalize_sum(wt0.ne, wt0.e);
    if (verbose) { wt_table_print(stderr, wname0, wt0.ne, wt0.e, 0); }
    
    uint32_t r1 = (uint32_t)ceil(sqrt(4*var1));
    uint32_t n1 = 2*r1 + 1;
    double_vec_t wt1 = wt_table_triangular_make(n1);
    char *wname1 = jsprintf("triangular(n=%d)", n1);
    if (verbose) { wt_table_print(stderr, wname1, wt1.ne, wt1.e, 0); }
    
    (*wt0_P) = wt0; (*wname0_P) = wname0;
    (*wt1_P) = wt1; (*wname1_P) = wname1;
  }
