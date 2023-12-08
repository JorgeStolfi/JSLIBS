#ifndef wt_table_binomial_H
#define wt_table_binomial_H

/* Weight tables with binomial profile. */
/* Last edited on 2023-11-25 12:05:48 by stolfi */

#define wt_table_binomial_H_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <argparser.h>

void wt_table_binomial_fill(int32_t n, double wt[], int32_t *stride_P);
  /* Stores into {wt[0..n-1]} a weight table with binomial profile.
    Namely, sets {wt[i]} to {choose(i,n-1)} for {i} in {0..n-1}.
    
    Note that the sum of all entries is {2^n}. The variance will be
    {(n-1)/4}, so the standard deviation will be {sqrt(n-1)/2)}. The
    table will be a partition of constant with stride 1, and (if {n>=2}) 
    also with stride 2. All these properties will hold exactly {n}
    up to {40} at least. */

double_vec_t wt_table_binomial_make(int32_t n);
  /* Allocates a new {double_vec_t} {wt} with {wt.ne = n} elements and fills {wt.e[0..n-1]} 
    with {wt_table_binomial_fill(n,wt.e,NULL)}. */

#endif
