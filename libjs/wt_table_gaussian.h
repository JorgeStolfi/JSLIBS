#ifndef wt_table_gaussian_H
#define wt_table_gaussian_H

/* Weight tables with Gaussain profile */
/* Last edited on 2023-11-25 12:09:15 by stolfi */

#define wt_table_gaussian_H_COPYRIGHT \
  "Copyright © 2023  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <argparser.h>

void wt_table_gaussian_fill(int32_t n, double sigma, double wt[], int32_t *stride_P);
  /* Stores into {wt[0..n-1]} a window weight table derived from a
    Gaussian distribution {G(z)} with mean {n/2} and standard deviation
    {sigma}.  The length {n} must be positive.
    
    More precisely, {wt[i]} will be the integral of {G(z)} over
    the interval {[i _ i+1]}.  The intergral in {(-oo _ 0]} and 
    {[n _ +oo)} is ignored.  
    
    In particular, if {sigma} is zero, the function {G} is assumed
    to be a Dirac impulse with integral 1. Then, if {n} is odd, {wt[n/2]} will be 1.0;
    if {n} is even, {wt[n/2-1]} and {wt[n/2]} will be 0.5; in either case
    all the other elements will be zero.
    
    If {sigma} is zero or {n} is 1 or 2, the table will be a partition
    of constant, with max stride 1 or 2. Otherwise the table is not
    exactly a partition of a constant (so the max stride is technically
    0) but is approximately so if {n} is large compared to {sigma} and
    replicated with a stride that is not too large compared to {sigma}.
    In any case, if {stride_P} not {NULL}, the max stride (0, 1, or 2)
    will be returned in {*stride_P}. */

double_vec_t wt_table_gaussian_make(int32_t n, double sigma, double maxLoss);
  /* Builds a weight table {wt} with odd length derived from a Gaussian
    distribution {G(z)} with standard deviation {sigma}.
    
    If {n} is not zero, the procedure makes {wt} a {double_vec_t} 
    with size {wt.ne = n}.
    
    If {n} is zero, the procedure chooses an odd table size {wt.ne =
    2*r+1} so that the total omitted (out-of-window) probability (the
    integral of {G(z)} outside the interval {[-(r+1/2) _ +(r+1/2)]})
    does not exceed {maxLoss}.
    
    In either case, the table is filled with 
    {wt_table_gaussian_fill(n,sigma,wt.e, NULL)}. */

/* AUXILIARY FUNCTIONS */

double wt_table_gaussian_loss(int32_t n, double sigma);
  /* Returns the weight that is lost if a table with {n} entries is
    used to represent a Gaussian weight distribution with standard
    deviation {sigma}. */

double wt_table_gaussian_entry(int32_t n, int32_t k, double sigma);
  /* Computes the entry with position {k} in a table with {n}
    elements, where each element is the integral of a truncated
    Gaussian weight distribution with standard deviation {sigma}
    in a unit sub-interval of the interval {[-n/2 _ +n/2]}. The
    position {k} ranges in {0..n-1}. The table is normalized
    to unit sum. 
    
    !!! Inefficient -- should reuse the shared {erf} calls !!! */

#endif
