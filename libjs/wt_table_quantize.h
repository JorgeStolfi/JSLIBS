#ifndef wt_table_quantize_H
#define wt_table_quantize_H

/* Quantiziing weight tables. */
/* Last edited on 2023-11-04 18:12:49 by stolfi */

#define wt_table_quantize_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

int32_t wt_table_quantize(int32_t n, double wd[], int32_t wi_min, int32_t wi_sum, int32_t wi[]);
  /* Stores into {wi[0..n-1]} an copy of the weights {wd[0..n-1]} scaled and rounded
    to integers. The input weights may be negative or zero. 
    
    If {wi_min} is nonzero, the procedure will ensure that {|wi[k]|}
    will be at least {wi_min} whenever {wd[k]} is nonzero. Otherwise
    some non-zero weights {wd[k]} may be rounded to zero. In any case,
    if {wi[k]} is not zero, it will have the same sign as {wd[k]}.
    
    If {wi_sum} is not 0, the scaling will be chosen so that the sum of
    the absolute values of {wi[0..n-1]} will be exactly {wi_sum}.
    
    If {wi_sum} is zero and {wi_min} is nonzero, the scale factor will
    be chosen so that the smallest non-zero entries of {wd[0..n-1]} are
    mapped to {wi_min}.
    
    However, if all input weights {wd[0..n-1]} are zero, the output
    weights {wi[0..n-1]} will be all zero, regardless of the parameters.
    
    Returns the sum of the absolute values of {wi_min[0..n-1]}.
    
    The value of {wi_sum} must be non-negative and at most
    {wt_table_quantize_WI_SUM_MAX}. The caller must make sure {wi_min}
    is not so large that the sum of the absolute values of {wi[0..n-1]}
    will exceed {wt_table_quantize_WI_SUM_MAX}. If both {wi_sum} and
    {wi_min} are zero, the procedure assumes {wi_sum =
    wt_table_quantize_WI_SUM_MAX}. */

#define wt_table_quantize_WI_SUM_MAX (INT32_MAX/2 -1)
  /* Maximum value of the {wi_sum} parameter of {wt_table_quantize}. */

#endif
