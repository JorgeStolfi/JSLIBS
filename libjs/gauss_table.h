#ifndef gauss_table_H
#define gauss_table_H

/* Tools for Gaussian kernel tables */
/* Last edited on 2024-11-19 05:59:17 by stolfi */

/* !!! Merge into {wt_table_gaussian} !!! */

#include <stdint.h>

#include <bool.h>
#include <gauss_bell.h>

#define gauss_table_HUGE_ARG (gauss_bell_HUGE_ARG)
  /* A value such that {|z| >= gauss_table_HUGE_ARG} implies 
    that the standard PDF at {z} underflows. */
     
#define gauss_table_BIG_ARG (gauss_bell_BIG_ARG)
  /* A value such that {|z| >= gauss_table_BIG_ARG} implies 
    that the standard PDF at {z} is less than {10^{-16}}. */
 
#define gauss_table_TINY_ARG (gauss_bell_TINY_ARG)
  /* A value such that {|z| <= gauss_table_TINY_ARG} implies
    that the standard PDF at {z} is greater than {1 - 10^{-16}}. */

double *gauss_table_make(uint32_t n, double avg, double dev, bool_t normSum, bool_t folded);
  /* Allocates a vector {w} with {n} elements and sets element {w[i]}
    to {g((i-avg)/dev)}, where {g(z)} is the Gaussian bell {exp(-z*z/2)}.
    
    If {folded} is true, the bell will be wrapped around the edges of
    the vector. That is, {w[i]} is actually set to {SUM {
    g((i-avg+k*n)/dev) : k in Z }}. The fold-over is considered not
    significant if the Gaussian is far from the table ends --- that is,
    if {BIG*dev <= avg <= n - BIG*dev} where {BIG = gauss_bell_BIG_ARG}.
    
    In particular, if {dev} is zero, entry {w[i]} is set to 1 if {i} is
    equal to {avg} (if {folded} is false) or congruent to {avg} modulo
    {n} (if {folded} is true), and to 0 otherwise. If {dev} is {+INF},
    all entries are set to 1.
    
    These raw values are normalized before the procedure returns. If
    {normSum} is TRUE, the elements {w[0..n-1]} are scaled so that
    their sum is 1.0. Otherwise, if {folded} is true, they are scaled 
    so that the weight at index {avg} (interpolated if {avg} is fractional,
    and including any wrap-around effects) is 1.0.  Note if {folded} is 
    false the value at {avg} is already 1.0. */

#define gauss_table_BIG_DEV (1.5)
  /* A value such that {dev >= h*gauss_table_BIG_DEV} implies
    that the Gaussian bell folded over {[0_h]} is practically constant. */

double gauss_table_folded_bell(double z, double dev, uint32_t n);
  /* Returns the sum of the Gaussian bell {exp(-(x/dev)^2/2)}
    evaluated at all arguments {x} that are congruent to {z} modulo {n}.
    
    As special cases, if {dev} is zero, returns 1 when {z} is
    congruent to 0 modulo {n}, and 0 otherwise. If {dev} is {+INF}, or
    greater than {n*gauss_table_BIG_DEV}, returns 1. If {n} is zero,
    assumes {n = +oo}, meaning no fold-over. Note that the result is
    not normalized, by any criterion. */

#endif
