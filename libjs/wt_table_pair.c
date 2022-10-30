/* See wt_table_pair.h */
/* Last edited on 2022-10-30 19:32:03 by stolfi */

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
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  )   
  { double maxLoss = 1.0e-3;
    double sigma0 = sqrt(var0);
    bool_t norm = TRUE;
    (*wtb0) = wt_table_make_gaussian(sigma0, maxLoss, norm);
    (*wname0) = NULL;
    asprintf(wname0, "gaussian(sigma = %8.6f)", sigma0);
    if (verbose) { wt_table_print(stderr, (*wname0), wtb0->ne, wtb0->e, 0); }
    
    double sigma1 = sqrt(var1);
    (*wtb1) = wt_table_make_gaussian(sigma1, maxLoss, norm);
    (*wname1) = NULL;
    asprintf(wname1, "gaussian(sigma = %8.6f)", sigma1);
    if (verbose) { wt_table_print(stderr, (*wname1), wtb1->ne, wtb1->e, 0); }
  }
   
void wt_table_pair_make_binomial
  ( double var0,
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  )   
  { bool_t norm = TRUE;
    
    /* Variance of the distribution {w[k] = choose(n,k)} is {n/4 = r/2}, so: */
    int r0 = (int)floor(2*var0 + 0.5);
    (*wtb0) = wt_table_make_binomial(r0, norm);
    (*wname0) = NULL;
    asprintf(wname0, "binomial(n = %d, sigma = %8.6f)", 2*r0, sqrt(0.5*r0));
    if (verbose) { wt_table_print(stderr, (*wname0), wtb0->ne, wtb0->e, 0); }
    
    int r1 = (int)floor(2*var1 + 0.5);
    (*wtb1) = wt_table_make_binomial(r1, norm);
    (*wname1) = NULL;
    asprintf(wname1, "binomial(n = %d, sigma = %8.6f)", 2*r1, sqrt(0.5*r1));
    if (verbose) { wt_table_print(stderr, (*wname1), wtb1->ne, wtb1->e, 0); }
  }
   
void wt_table_pair_make_triangular
  ( double var0,
    double_vec_t *wtb0, 
    char **wname0, 
    double var1,
    double_vec_t *wtb1, 
    char **wname1,
    bool_t verbose
  )   
  { bool_t norm = TRUE;
    
    /* Variance of the distribution {w[k] = r + 1 - |r-k|} is {~r^2/4}, so: */
    int r0 = (int)floor(sqrt(4*var0) + 0.5);
    (*wtb0) = wt_table_make_triangular(r0, norm);
    (*wname0) = NULL;
    asprintf(wname0, "triangular(r = %d)", r0);
    if (verbose) { wt_table_print(stderr, (*wname0), wtb0->ne, wtb0->e, 0); }
    
    int r1 = (int)floor(sqrt(4*var1) + 0.5);
    (*wtb1) = wt_table_make_triangular(r1, norm);
    (*wname1) = NULL;
    asprintf(wname1, "triangular(r = %d)", r1);
    if (verbose) { wt_table_print(stderr, (*wname1), wtb1->ne, wtb1->e, 0); }
  }
