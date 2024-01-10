#ifndef gauss_table_H
#define gauss_table_H

/* Tools for Gaussian kernel tables */
/* Last edited on 2011-05-09 13:22:14 by stolfi */

#include <bool.h>

#define gauss_table_BIG_ARG (8.6)
  /* A value such that {z >= gauss_table_BIG_ARG} implies
    {exp(-z^2/2) <= 1.0e-16}. */
 
#define gauss_table_TINY_ARG (1.4e-8)
  /* A value such that {0 < z <= gauss_table_TINY_ARG} implies
    {exp(-z^2/2) >= 1 - 1.0e-16}. */
 
#define gauss_table_BIG_DEV (1.5)
  /* A value such that {dev >= n*gauss_table_BIG_DEV} implies
    that the Gaussian bell folded over {[0_n]} is constant. */

double *gauss_table_make(int n, double avg, double dev, bool_t normSum);
  /* Allocates a vector {w} with {n} elements and sets element {w[i]}
    to {g((i-avg)/dev)}, where {g(z)} is the Gaussian bell {exp(-z*z/2)}.
    
    The bell will be wrapped around the edges of the vector. That is,
    {w[i]} is actually set to {SUM { g((i-avg+k*n)/dev) : k in Z }}.
    The fold-over is not significant if the Gaussian is far from the
    table ends --- that is, if {BIG*dev <= avg <= n - BIG*dev} where
    {BIG = gauss_table_BIG_ARG}.
    
    In particular, if {dev} is zero, entry {w[i]} is set to 1 if {i}
    is congruent to {avg} modulo {n}, and to 0 otherwise. If {dev} is
    {+INF}, all entries are set to 1.
    
    These raw values are normalized before the procedure returns. If
    {normSum} is TRUE, the elements {w[0..n-1]} are scaled so that
    their sum is 1.0. Otherwise they are scaled so that the weight at
    index {avg} (interpolated if {avg} is fractional, and including
    any wrap-around effects) is 1.0. */

double gauss_table_folded_bell(double z, double dev, int n);
  /* Returns the sum of the Gaussian bell {exp(-(x/dev)^2/2)} evaluated 
    at all arguments {x} that are congruent to {z} modulo {n}.
    
    As special cases, if {dev} is zero, returns 1 when {z} is
    congruent to 0 modulo {n}, and 0 otherwise. If {dev} is {+INF}, or
    greater than {n*gauss_table_BIG_DEV}, returns 1. If {n} is zero,
    assumes {n = +oo}, meaning no fold-over. Note that the result is
    not normalized, by any criterion. */

double gauss_table_bell(double z, double dev);
  /* The unnormalized Gaussian bell {exp(-(x/dev)^2/2)}, but safer/faster 
     for big and small args. */

#endif
