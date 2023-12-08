#ifndef wt_table_uniform_H
#define wt_table_uniform_H

/* Uniform (constant) weight tables. */
/* Last edited on 2023-11-25 12:04:13 by stolfi */

#define wt_table_uniform_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <argparser.h>

void wt_table_uniform_fill(int32_t n, double val, double wt[], int32_t *stride_P);
  /* Fills {wt[0..n-1]} with the value {val}.  The table size {n} and the 
    value {val} must be positive.
    
    The table is a partition of constant (PoC) when used with stride {n}
    or any divisor thereof. */

double_vec_t wt_table_uniform_make(int32_t n, double val);
  /* Allocates a new {double_vec_t} {wt} with {wt.ne = n} elements and fills {wt.e[0..n-1]} 
    with {wt_table_uniform_fill(n, val, wt.e, NULL)}. */

#endif
