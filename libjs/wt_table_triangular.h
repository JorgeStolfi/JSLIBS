#ifndef wt_table_triangular_H
#define wt_table_triangular_H

/* Triangular weight tables*/
/* Last edited on 2023-11-25 12:06:05 by stolfi */

#define wt_table_triangular_H_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <argparser.h>

void wt_table_triangular_fill(int32_t n, double wt[], int32_t *stride_P);
  /* Stores into {wt[0..n-1]} a window weight table with triangular profile.  
  
    More precisely, sets {wt[i] = i + 1} for {i} in {0..(n-1)/2}, and
    {wt[i] = wt[n-1-i]} for the other elements. Thus, if {n} is odd, the
    central element {wt[r]} will be set to {(n+1)/2}, and the table will
    be the convolution of two uniform distributions on {r+1} elements.
    If {n} is even the two centermost elements {w[n/2-1]} and {w[n/2]}
    will be set to {n/2}.
    
    In any case, the values will be symmetric; the elements {wt[-1]} and
    {wt[n]}, if they existed, would be zero, and {wt[0]} and {wt[n-1]}
    will be 1.
    
    The table will be a partition of constant (PoC) with max stride
    {(n+1)/2}. If {stride_P} is not {NULL}, the max stride will be
    returned in {*stride_P}.    
    
    !!! What is the variance? !!! */

double_vec_t wt_table_triangular_make(int32_t n);
  /* Allocates a new {double_vec_t} {wt} with {wt.ne = n} elements and
    fills {wt.e[0..n-1]} with {wt_table_triangular_fill(n,wt.e,NULL)}. */

#endif
