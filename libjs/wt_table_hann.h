#ifndef wt_table_hann_H
#define wt_table_hann_H

/* Weight table with Hann (shifted cosine) profile and optional flat center. */
/* Last edited on 2023-11-25 12:04:35 by stolfi */

#define wt_table_hann_H_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <argparser.h>

void wt_table_hann_fill(int32_t n, double flat, double wt[], int32_t *stride_P);
  /* Stores into {wt[0..n-1]} a window weight table with a 
    Hann-like sinusoidal rise to 1.0, a flat (constant 1.0) section,
    and sinusoidal drop to 0.  The length {n} must be positive.
    
    More precisely, sets {wt[i] = 1.0} if {r <= i <= n-1-r}, to
    {0.5*(1+cos(pi*(i-r)/(r+1)))} if {0 <= i < r}, and to {wc[n-1-i]} if
    {i > n-1-r}; where {r = 0.5*(1-flat)*(n-1)}. 
    
    In particular, if {flat} is zero (so that {r = 0.5*(n-1)}, the
    result is the discrete Hann (shifted cosine) distribution.  
    If {flat} is 1 (so that {r = 0}), the table uniform equal to 1.
    
    The table is a partition of constant (PoC) and only if {r} is an
    integer (in particular, if {flat} is 1, or {n} is odd and {flat} is
    0), in which case its max PoC stride will be {n-r}. Otherwise the
    max PoC stride is 0 by convention. If {stride_P} is not {NULL},
    {*stride_P} is set to the max PoC stride. */

double_vec_t wt_table_hann_make(int32_t n, double flat);
  /* Allocates a new {double_vec_t} {wt} with {wt.ne = n} elements and fills {wt.e[0..n-1]} 
    with {wt_table_hann_fill(n,flat,wt.e,NULL)}. */

#endif
