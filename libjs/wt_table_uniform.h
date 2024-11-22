#ifndef wt_table_uniform_H
#define wt_table_uniform_H

/* Uniform (constant) weight tables. */
/* Last edited on 2024-11-16 10:31:38 by stolfi */

#define wt_table_uniform_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <argparser.h>

void wt_table_uniform_fill(uint32_t n, double val, double wt[], uint32_t *stride_P);
  /* Fills {wt[0..n-1]} with the value {val}.  The table size {n} and the 
    value {val} must be positive.
    
    The table is a partition of constant (PoC) when used with stride {n}
    or any divisor thereof. */

double_vec_t wt_table_uniform_make(uint32_t n, double val);
  /* Allocates a new {double_vec_t} {wt} with {wt.ne = n} elements and fills {wt.e[0..n-1]} 
    with {wt_table_uniform_fill(n, val, wt.e, NULL)}. */

#endif
