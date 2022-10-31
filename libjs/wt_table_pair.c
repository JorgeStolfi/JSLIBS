/* See wt_table_pair.h */
/* Last edited on 2022-10-31 14:11:29 by stolfi */

#define wt_table_pair_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <vec.h>
#include <bool.h>
#include <affirm.h>

#include <wt_table.h>
#include <wt_table_pair.h>

/* INTERNAL PROTOTYPES */
   
/* IMPLEMENTATIONS */

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
    bool_t norm = TRUE;
    double_vec_t wt0 = wt_table_make_gaussian(sigma0, maxLoss, norm);
    char *wname0 = NULL;
    int32_t n0 = wt0.ne;
    asprintf(&wname0, "gaussian(n=%d,sigma=%8.6f)", n0, sigma0);
    if (verbose) { wt_table_print(stderr, wname0, wt0.ne, wt0.e, 0); }
    
    double sigma1 = sqrt(var1);
    double_vec_t wt1 = wt_table_make_gaussian(sigma1, maxLoss, norm);
    char *wname1 = NULL;
    int32_t n1 = wt1.ne;
    asprintf(&wname1, "gaussian(n=%d,sigma=%8.6f)", n1, sigma1);
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
  { bool_t norm = TRUE;
    
    /* Variance of the distribution {w[k] = choose(n,k)} is {n/4 = r/2}, so: */
    int32_t r0 = (int32_t)ceil(2*var0);
    int32_t n0 = 2*r0 + 1;
    double_vec_t wt0 = wt_table_make_binomial(n0, norm);
    char *wname0 = NULL;
    asprintf(&wname0, "binomial(n=%d)", n0);
    if (verbose) { wt_table_print(stderr, wname0, wt0.ne, wt0.e, 0); }
    
    int32_t r1 = (int32_t)ceil(2*var1);
    int32_t n1 = 2*r1 + 1;
    double_vec_t wt1 = wt_table_make_binomial(n1, norm);
    char *wname1 = NULL;
    asprintf(&wname1, "binomial(n=%d)", n1);
    if (verbose) { wt_table_print(stderr, wname1, wt1.ne, wt1.e, 0); }
    
    (*wt0_P) = wt0; (*wname0_P) = wname0;
    (*wt1_P) = wt1; (*wname1_P) = wname1;
  }
   
void wt_table_pair_make_triangular
  ( double var0,
    double_vec_t *wt0, 
    char **wname0_P, 
    double var1_P,
    double_vec_t *wt1_P, 
    char **wname1_P,
    bool_t verbose
  )   
  { bool_t norm = TRUE;
    
    /* Variance of the distribution {w[k] = r + 1 - |r-k|} is {~r^2/4}, so: */
    int32_t r0 = (int32_t)ceil(sqrt(4*var0));
    int32_t n0 = 2*r0 + 1;
    (*double_vec_t wt0 = wt_table_make_triangular(n0, norm);
    char *wname0 = NULL;
    asprintf(&wname0, "triangular(n=%d)", n0);
    if (verbose) { wt_table_print(stderr, wname0, wt0.ne, wt0.e, 0); }
    
    int32_t r1 = (int32_t)ceil(sqrt(4*var1));
    int32_t n1 = 2*r1 + 1;
    double_vec_t wt1 = wt_table_make_triangular(n1, norm);
    char *wname1 = NULL;
    asprintf(&wname1, "triangular(n=%d)", n1);
    if (verbose) { wt_table_print(stderr, wname1, wt1.ne, wt1.e, 0); }
    
    (*wt0_P) = wt0; (*wname0_P) = wname0;
    (*wt1_P) = wt1; (*wname1_P) = wname1;
  }
