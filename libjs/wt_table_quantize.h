#ifndef wt_table_quantize_H
#define wt_table_quantize_H

/* Quantiziing weight tables. */
/* Last edited on 2024-11-19 09:08:24 by stolfi */

#define wt_table_quantize_H_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

uint64_t wt_table_quantize
  ( uint32_t n,
    double wf[],
    uint64_t wia_sum,
    bool_t keep_nz,
    int32_t wi[]
  );
  /* Stores into {wi[0..n-1]} a copy of the weights {wf[0..n-1]}
    scaled and rounded to integers, so that the sum of their absolute 
    values is about {wia_sum}.  The input weights may be negative or zero. 
    
    If {keep_nz} is true,the procedure will make sure that each output
    weight {wi[i]} is nonzero whenever {wf[i]} is nonzero. If {keep_nz}
    is false, some nonzero input wieights {wf[i]} may be rounded to zero.
    
    More precisely, the procedure initially sets a scale factor {s} to
    {wia_sum/wfa_sum}, where {wfa_sum} is the sum of the absolute values of
    all input weights {wf[0..n-1]}.
    
    Then it sets each {wi[i]} to {s*wf[i]} rounded to the nearest integer.
    Then, if {keep_nz} is true, whenever {wi[i]=0} and {wf[i]!=0}, {wi[i]} 
    is set to the sign ({-1} or {+1}) of {wf[i]}.  
    
    Then {s} is adjusted by binary search, and the scaling and rounding
    of the previous paragraph is repeated, until the sum of the absolute values
    of the {wi[i]} is as close as possible to {wia_sum}.
    
    As a special case, if the input weights are all zero, the outputs
    {wi[0..n-1]} will be all zeros too, and {wia_sum} will be ignored.
    
    Also, the outputs {wi[0..n-1]} will be all zeros if {wia_sum} 
    is zero or too small and {keep_nz} is false.
    
    The procedure returns the actual sum of the absolute values of the
    rounded weights {wi[0..n-1]}.
    
    The parameter {wia_sum} must not exceed {wt_table_quantize_WIA_SUM_MAX},
    and it will be reduced as needed to ensure that no output weight
    exceeds {wt_table_quantize_WIA_MAX} in absolute value */
 
#define wt_table_quantize_WIA_SUM_MAX (1LU << 50)
  /* Maximum abs value of the {wia_sum} parameter of {wt_table_quantize}. */
    
#define wt_table_quantize_WIA_MAX (UINT32_MAX/8)
  /* Maximum absolute value allowed for the output weights. */
  
#define wt_table_quantize_N_MAX (1U << 26)
  /* Maximum value of {n} for {wt_table_quantize}. */

#endif
